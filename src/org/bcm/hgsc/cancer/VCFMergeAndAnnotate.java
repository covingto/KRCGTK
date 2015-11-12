package org.bcm.hgsc.cancer;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.lang3.ArrayUtils;
import org.bcm.hgsc.cancer.MultiVCFReader.AlleleContainer;
import org.bcm.hgsc.cancer.utils.CARNACSampleGenotyper;
import org.bcm.hgsc.cancer.utils.SampleGenotyper;
import org.bcm.hgsc.utils.AlleleResolver;
import org.bcm.hgsc.utils.AlleleResolver.AlleleSet;
import org.bcm.hgsc.utils.BAMInterface;
import org.bcm.hgsc.utils.BAMUtils;
import org.bcm.hgsc.utils.BAMUtils.ConformedRead;
import org.bcm.hgsc.utils.Settings;
import org.bcm.hgsc.utils.Utils;

//import org.broadinstitute.variant.variantcontext.writer.

/**
 * Merge and annotate VCF files. This merges multiple VCF files together for
 * both tumor and normal samples, extracts relevant sequence information from
 * the BAM files and writes a new VCF file with the information. This method
 * extracts all required information for making a call in the CARNAC variant
 * filter.
 * 
 * @author covingto
 * 
 */
public class VCFMergeAndAnnotate {
	private static Logger log = Logger.getLogger("");
	public static List<ThreadedAlleleResolver> workers = null;

	private static enum State {
		READING, PROCESSING, WRITING, DONE
	};

	private class StateMonitor {
		private State state = State.READING;

		StateMonitor() {
		}

		public void setProcessing() {
			this.state = State.PROCESSING;
		}

		public void setWriting() {
			this.state = State.WRITING;
		}

		public void setDone() {
			this.state = State.DONE;
		}

		public State getState() {
			return this.state;
		}
	}

	private class ThreadedAlleleContainerLoader implements Runnable {
		private final MultiVCFReader reader;
		private final BlockingQueue<AlleleContainer> container;
		private final BlockingQueue<AlleleContainer> writerContainer;
		private boolean done = false;
		private final StateMonitor monitor;

		ThreadedAlleleContainerLoader(MultiVCFReader reader,
				BlockingQueue<AlleleContainer> container,
				BlockingQueue<AlleleContainer> writerContainer,
				StateMonitor monitor) {
			this.reader = reader;
			this.container = container;
			this.writerContainer = writerContainer;
			this.monitor = monitor;
		}

		@Override
		public void run() {
			while (this.reader.hasNext()) {
				try {
					AlleleContainer ac = this.reader.nextAlleleSet();
					synchronized (log) {
						log.log(Level.FINE,
								"Adding allele container " + ac.toString());
					}
					container.put(ac);
					writerContainer.put(ac);
				} catch (InterruptedException e) {
					if (this.reader.hasNext()) {
						log.log(Level.SEVERE,
								"Reader interrupted while adding allele to pool.");
					}
					break;
				}
			}
			this.done = true;
			synchronized (log) {
				log.log(Level.INFO, "Done loading allele containers");
			}
			return;
		}

		public boolean done() {
			return done;
		}
	}

	private class ThreadedAlleleWriter implements Runnable {
		private final VariantContextWriter writer;
		private final BlockingQueue<AlleleContainer> alleleContainer; // this is the writer container
		private final Map<AlleleContainer, VariantContext> resultMap;
		private boolean done = false;
		private final StateMonitor monitor;

		public ThreadedAlleleWriter(VariantContextWriter writer,
				BlockingQueue<AlleleContainer> alleleContainer,
				Map<AlleleContainer, VariantContext> resultMap,
				StateMonitor monitor) {
			this.writer = writer;
			this.alleleContainer = alleleContainer;
			this.resultMap = resultMap;
			this.monitor = monitor;
		}

		@Override
		public void run() {
			while (!Thread.currentThread().isInterrupted()
					&& !monitor.getState().equals(State.WRITING)) {
				try {
					/**
					 * @note it is important that the writer does poll instead of take
					 * because there is the possibility that the alleleContainer is empty before 
					 * @ref State.WRITING can be set.
					 */
					final AlleleContainer ac = this.alleleContainer.poll(30, TimeUnit.SECONDS);
					if (ac != null){
						resolve(ac);
					}
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
			}
			final List<AlleleContainer> finalContainers = new ArrayList<AlleleContainer>();
			this.alleleContainer.drainTo(finalContainers);
			for (final AlleleContainer ac : finalContainers) {
				resolve(ac);
			}
			this.done = true;
		}
		
		private void resolve(AlleleContainer ac) {
			log.fine("Resolving allele container " + ac.toString());
			while (!resultMap.containsKey(ac) && !monitor.getState().equals(State.WRITING)) {
				try {
					log.finer("Waiting for allele " + ac.toString() + " to become available");
					Thread.sleep(1000);
					
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
			}
			final VariantContext vc = resultMap.get(ac);
			if (vc != null) {
				log.fine("Adding variant to writer");
				writer.add(vc);
				log.fine("Variant added to writer");
			} else {
				log.log(Level.WARNING, "Failed to find variant context " + ac.toString());
			}
			synchronized (resultMap) {
				resultMap.remove(ac);
			}
		}

		public boolean done() {
			return done;
		}
	}

	private class SingleThrowThreadedAlleleResolver implements Runnable {
		private final int includeFlag;
		private final int excludeFlag;
		private final List<BAMInterface> baminterfaces;
		// private final BlockingQueue<VariantContext> variantContextQueue; //
		// will be used only once to add the variant
		private final AlleleResolver.ResolutionType resolution;
		private final SampleGenotyper genotyper;
		private final File fastafile;
		private final int padding;
		// private final VariantContextWriter writer;
		private final AlleleContainer alleleContainer;
		private final Map<AlleleContainer, VariantContext> resultMap;
		private final StateMonitor monitor;

		/**
		 * Class initialization should contain all information required to run
		 * the method. Running is just a pass through
		 * 
		 * @param monitor
		 */
		public SingleThrowThreadedAlleleResolver(int padding,
				AlleleResolver.ResolutionType resolution,
				SampleGenotyper genotyper, List<BAMInterface> baminterfaces,
				Map<AlleleContainer, VariantContext> resultMap, File fastafile,
				AlleleContainer container, StateMonitor monitor,
				int f, int F) {
			// this.alleleContainer = alleleContainer;
			this.resolution = resolution;
			this.genotyper = genotyper;
			this.baminterfaces = baminterfaces;
			// this.variantContextQueue = variantContextQueue;
			this.padding = padding;
			this.fastafile = fastafile;
			this.resultMap = resultMap;
			this.alleleContainer = container;
			this.monitor = monitor;
			this.includeFlag = f;
			this.excludeFlag = F;
		}

		@Override
		public void run() {
			log.log(Level.FINEST, "Allele processing");
			try {
				final VariantContext newv = processAlleleContainer(alleleContainer);
				synchronized (this.resultMap) {
					log.log(Level.FINEST, "Adding allele result to resultMap");
					this.resultMap.put(alleleContainer, newv);
					log.log(Level.FINEST, "Allele added to resultMap");
				}
			} catch (Exception e) {
				log.log(Level.SEVERE, "Caught exception while processing allele", e);
				synchronized (this.resultMap) {
					log.log(Level.FINEST, "Adding null result to resultMap");
					this.resultMap.put(alleleContainer, null);
					log.log(Level.FINEST, "Allele null added to resultMap");
				}
			}
		}

		private VariantContext processAlleleContainer(
				AlleleContainer alleleContainer) {
			VariantContext newv = null;
			final int start = alleleContainer.getStart() - padding;
			final int end = alleleContainer.getEnd() + padding;
			try {
				// log.log(Level.INFO, "Building gene holders");
				final List<ConformedRead> allReads = new ArrayList<ConformedRead>(
						1000);
				final Map<String, List<ConformedRead>> sampleReads = new HashMap<String, List<ConformedRead>>();
				AlleleSet alleles = null;
				try {
					// IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(this.fastafile);
					for (BAMInterface bi : baminterfaces) {
						final SamReader sam = bi.getSamfilereader();
						
						final List<ConformedRead> reads = new ArrayList<ConformedRead>(1000);
						sampleReads.put(bi.getSampleName(), reads);
						Iterator<ConformedRead> cri = BAMUtils.getConformedReadsIterator(
								sam, alleleContainer.getChr(),
								start, end,
								this.includeFlag, this.excludeFlag,
								new IndexedFastaSequenceFile(this.fastafile));
						while (cri.hasNext()){
							final ConformedRead cr = cri.next();
							// if null, we don't want it
							// if not covering slice, we don't want it
							// if covering undefined (dot) bases (like introns), we don't want it
							if (cr == null || cr.readStart() > start || cr.readEnd() < end) {
								continue;
							} 
							if (ArrayUtils.contains(cr.getReadAtGenomicRange(start, end).bytes, BAMUtils.dot)){
								// log.log(Level.FINE, "Read contains 'dot' in " + start + "-" + end + ":\n" + cr.toString());
								continue;
							}
							reads.add(cr);
						}
						
						allReads.addAll(reads);
						sam.close();
					}
					// synchronized (log) {
					// log.log(Level.FINE, "Simplifying allele");
					// }
					alleles = AlleleResolver.resolveAlleles(allReads,
							alleleContainer.getChr(),
							alleleContainer.getStart() - padding,
							alleleContainer.getEnd() + padding, resolution,
							new IndexedFastaSequenceFile(this.fastafile),
							this.padding).simplify();

					// add the info to the variant context builder
					if (alleles == null || alleles.getAlleles().size() == 0) {
						// in case the alleles are null, we still want to output
						// something, so we will default to the standard set
						synchronized (log) {
							log.log(Level.FINE,
									"Adding raw contexts with no genotypes since no alleles returned from allele resolver.");
						}
						for (VariantContext vc : alleleContainer
								.getVariantContexts()) {
							VariantContextBuilder vcbuilder = new VariantContextBuilder();
							vcbuilder.alleles(vc.getAlleles());
							synchronized (log) {
								log.log(Level.FINE,
										"Adding null context to writer");
							}
							newv = vcbuilder.make();
						}

					} else {
						// build the attribute map
						Map<String, Object> varAttributes = new TreeMap<String, Object>();
						varAttributes
								.put("OC", alleleContainer.getCallString());
						VariantContextBuilder vcbuilder = new VariantContextBuilder();

						vcbuilder.alleles(alleles.getAlleles());
						List<Genotype> genotypes = new ArrayList<Genotype>(
								sampleReads.size());
						// do the genotyping
						for (String k : sampleReads.keySet()) {
							final List<ConformedRead> reads = sampleReads
									.get(k);
							genotypes
									.add(genotyper.genotype(k, alleles, reads));
						}
						// TODO: filter the genotypes for quality
						vcbuilder
								.loc(alleleContainer.getChr(),
										alleles.getStart(), alleles.getEnd())
								.attributes(varAttributes).genotypes(genotypes);
						newv = vcbuilder.make();
					}
				} catch (Exception e) {
					if (alleles != null){
						log.log(Level.SEVERE, "Hit error in determining genotype: "
							+ alleleContainer.toString() + "\n" + alleles.toString(), e);
					} else {
						log.log(Level.SEVERE, "Hit error in determining genotype: "
								+ alleleContainer.toString() + "\n" + "Alleles were NULL", e);
					}
				} 
			} catch (Exception e) {
				log.log(Level.SEVERE,
						"Error in resolving allele, will add null to the allele map", e);
			}
			log.log(Level.FINEST, "Allele processed");
			return newv;
		}
	}
	
	private class ThreadedAlleleResolver implements Runnable {
		private final int includeFlag;
		private final int excludeFlag;
		private final List<BAMInterface> baminterfaces;
		// private final BlockingQueue<VariantContext> variantContextQueue; //
		// will be used only once to add the variant
		private final AlleleResolver.ResolutionType resolution;
		private final SampleGenotyper genotyper;
		private final File fastafile;
		private final int padding;
		// private final VariantContextWriter writer;
		private final BlockingQueue<AlleleContainer> container;
		private final Map<AlleleContainer, VariantContext> resultMap;
		private final StateMonitor monitor;

		/**
		 * Class initialization should contain all information required to run
		 * the method. Running is just a pass through
		 * 
		 * @param monitor
		 */
		public ThreadedAlleleResolver(int padding,
				AlleleResolver.ResolutionType resolution,
				SampleGenotyper genotyper, List<BAMInterface> baminterfaces,
				Map<AlleleContainer, VariantContext> resultMap, File fastafile,
				BlockingQueue<AlleleContainer> container, StateMonitor monitor,
				int f, int F) {
			// this.alleleContainer = alleleContainer;
			this.resolution = resolution;
			this.genotyper = genotyper;
			this.baminterfaces = baminterfaces;
			// this.variantContextQueue = variantContextQueue;
			this.padding = padding;
			this.fastafile = fastafile;
			this.resultMap = resultMap;
			this.container = container;
			this.monitor = monitor;
			this.includeFlag = f;
			this.excludeFlag = F;
		}

		@Override
		public void run() {
			log.log(Level.FINE, "Starting allele processing");
			while (!Thread.currentThread().isInterrupted()
					&& !monitor.getState().equals(State.PROCESSING)
					&& !monitor.getState().equals(State.WRITING)) {
				// polling will empty the container unless it is already empty,
				// in which case we just keep going.  The loop is broken when the state monitor enters the processing phase
				// at which point the stack is drained.
				try {
					final AlleleContainer alleleContainer = this.container.take(); 
					final VariantContext newv = processAlleleContainer(alleleContainer);
					synchronized (this.resultMap) {
						this.resultMap.put(alleleContainer, newv);
					}
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
				
			}
			// drain the pool and resolve
			log.log(Level.INFO, "Draining the pool and processing.");
			final List<AlleleContainer> finalContainers = new ArrayList<AlleleContainer>();
			synchronized (this.container) {
				this.container.drainTo(finalContainers);
			}
			for (final AlleleContainer alleleContainer : finalContainers) {
				final VariantContext newv = processAlleleContainer(alleleContainer);
				synchronized (this.resultMap) {
					this.resultMap.put(alleleContainer, newv);
				}
			}
		}

		private VariantContext processAlleleContainer(
				AlleleContainer alleleContainer) {
			VariantContext newv = null;
			final int start = alleleContainer.getStart() - padding;
			final int end = alleleContainer.getEnd() + padding;
			try {
				// log.log(Level.INFO, "Building gene holders");
				final List<ConformedRead> allReads = new ArrayList<ConformedRead>(
						1000);
				final Map<String, List<ConformedRead>> sampleReads = new HashMap<String, List<ConformedRead>>();
				AlleleSet alleles = null;
				try {
					// IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(this.fastafile);
					for (BAMInterface bi : baminterfaces) {
						final SamReader sam = bi.getSamfilereader();
						final List<ConformedRead> reads = new ArrayList<ConformedRead>(1000);
						sampleReads.put(bi.getSampleName(), reads);
						Iterator<ConformedRead> cri = BAMUtils.getConformedReadsIterator(
								sam, alleleContainer.getChr(),
								start, end,
								this.includeFlag, this.excludeFlag,
								new IndexedFastaSequenceFile(this.fastafile));
						while (cri.hasNext()){
							final ConformedRead cr = cri.next();
							// if null, we don't want it
							// if not covering slice, we don't want it
							// if covering undefined (dot) bases (like introns), we don't want it
							if (cr == null || cr.readStart() > start || cr.readEnd() < end) {
								continue;
							} 
							if (ArrayUtils.contains(cr.getReadAtGenomicRange(start, end).bytes, BAMUtils.dot)){
								log.log(Level.FINE, "Read contains 'dot' in " + start + "-" + end + ":\n" + cr.toString());
								continue;
							}
							reads.add(cr);
						}
						
						allReads.addAll(reads);
						sam.close();
					}
					// synchronized (log) {
					// log.log(Level.FINE, "Simplifying allele");
					// }
					alleles = AlleleResolver.resolveAlleles(allReads,
							alleleContainer.getChr(),
							alleleContainer.getStart() - padding,
							alleleContainer.getEnd() + padding, resolution,
							new IndexedFastaSequenceFile(this.fastafile),
							this.padding).simplify();

					// add the info to the variant context builder
					if (alleles == null || alleles.getAlleles().size() == 0) {
						// in case the alleles are null, we still want to output
						// something, so we will default to the standard set
						synchronized (log) {
							log.log(Level.FINE,
									"Adding raw contexts with no genotypes since no alleles returned from allele resolver.");
						}
						for (VariantContext vc : alleleContainer
								.getVariantContexts()) {
							VariantContextBuilder vcbuilder = new VariantContextBuilder();
							vcbuilder.alleles(vc.getAlleles());
							synchronized (log) {
								log.log(Level.FINE,
										"Adding null context to writer");
							}
							newv = vcbuilder.make();
						}

					} else {
						// build the attribute map
						Map<String, Object> varAttributes = new TreeMap<String, Object>();
						varAttributes
								.put("OC", alleleContainer.getCallString());
						VariantContextBuilder vcbuilder = new VariantContextBuilder();

						vcbuilder.alleles(alleles.getAlleles());
						List<Genotype> genotypes = new ArrayList<Genotype>(
								sampleReads.size());
						// do the genotyping
						for (String k : sampleReads.keySet()) {
							final List<ConformedRead> reads = sampleReads
									.get(k);
							genotypes
									.add(genotyper.genotype(k, alleles, reads));
						}
						// TODO: filter the genotypes for quality
						vcbuilder
								.loc(alleleContainer.getChr(),
										alleles.getStart(), alleles.getEnd())
								.attributes(varAttributes).genotypes(genotypes);
						newv = vcbuilder.make();
					}
				} catch (Exception e) {
					if (alleles != null){
						log.log(Level.SEVERE, "Hit error in determining genotype: "
							+ alleleContainer.toString() + "\n" + alleles.toString(), e);
					} else {
						log.log(Level.SEVERE, "Hit error in determining genotype: "
								+ alleleContainer.toString() + "\n" + "Alleles were NULL", e);
					}
				} finally {
					// explicitly set items to null this might help the garbage
					// collector, anyway these things are in the young
					// generation anyway...
					// log.log(Level.INFO, "Clearing reads");
					allReads.clear();
					for (List<ConformedRead> v : sampleReads.values()) {
						v.clear();
					}
					sampleReads.clear();
				}
				// log.log(Level.FINE, "Variant resolution complete");
			} catch (Exception e) {
				log.log(Level.SEVERE,
						"Error in resolving allele, will add null to the allele map", e);
			}
			return newv;
		}
	}

	public static void main(String[] args) throws Exception {

		// build a VCF reader object, this returns VCF records in batches for processing through the site lookup,
		// VCF records are returned in batches separated by an indicated buffer size between the returned batch and the next batch available.
		// These batches are then converted to genotypes and each genotype is tested.

		Options options = new Options();
		Parser parser = new BasicParser();
		Option vcfFilesOption = new Option("v", "vcf", true, "VCF files from callers or other input (required)");
		vcfFilesOption.setArgs(Option.UNLIMITED_VALUES);
		options.addOption(vcfFilesOption);
		Option bamFileOption = new Option("bam", true, "BAM file in which to report alleles (required)");
		bamFileOption.setArgs(Option.UNLIMITED_VALUES);
		options.addOption(bamFileOption);
		Option sampleNameOption = new Option("samplename", true, "Sample names for BAM files, will lead to swapped data if these are not in order, must be same length as BAM files.");
		sampleNameOption.setArgs(Option.UNLIMITED_VALUES);
		options.addOption(sampleNameOption);
		//options.addOption("tumorBAM", true, "tumor BAM file");
		//options.addOption("normalBAM", true, "normal BAM file");
		//options.addOption("tumorName", true, "name of the tumor sample. default is to extract from BAM name [basename <tumorBAM> .bam]");
		//options.addOption("normalName", true, "name of the normal sample. default is to extract from BAM name [basename <normalBAM> .bam]");
		options.addOption("r", true, "indexed reference file (required)");
		options.addOption("b", true, "buffer for VCF allele lookup.  variants within the buffer will be considered jointly for genotyping and allele resolution. [20]");
		options.addOption("p", true, "padding to be applied around allele sets (based on buffer).  padding prevents expansion of alleles but does not prevent re-genotyping (nothing does)");
		options.addOption("minExp", false, "provide minimal expansion of alleles, alleles only expand when SNP is beside INDEL");
		options.addOption("nThreads", true, "number of threads in which to process variants");
		options.addOption("maxSize", true, "maximum size of event to report, deletions larger than this size will be omitted from the output [20]");
		options.addOption("o", true, "output vcf file.  this file will be a vcf formatted file, please add .vcf at the end of the file name [<tumorName>_<normalName>.vcf]");
		options.addOption("h", false, "print help");
		options.addOption("d", false, "Should the program be run in debug mode (outputs all)");
		options.addOption("f", true, "required flag (same as in samtools) defaults to 0");
		options.addOption("F", true, "filtering flag (same as in samtools) defaults to 1284 [unmapped, not primary alignment, duplicates]");
		options.addOption("minAlleleCount", true, "Do not report alleles below this minimum value.  " + 
				"In general it is not recomended to use this because the user should be filtering at a later step." + 
				"  However, for some technologies, the error mode is so high that a reasonable filter can be implemented.");
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line = parser.parse(options, args);
		Settings.debug = line.hasOption("d");

		ConsoleHandler handler = new ConsoleHandler();
		handler.setFormatter(Settings.defautlFormatter());
		log.addHandler(handler);
		if (Settings.debug){
			handler.setLevel(Level.ALL);
			log.setLevel(Level.ALL);
		} else {
			handler.setLevel(Level.INFO);
			log.setLevel(Level.INFO);
		}

		if (! line.hasOption("vcf") || ! line.hasOption("bam") || ! line.hasOption("r") || line.hasOption("h")){
			log.log( Level.SEVERE, "Required args not supplied.");
			formatter.printHelp("Wheeljack.jar", options);
			System.exit(10);
		}
		File[] 		vcffiles = Utils.filesFromStrings(line.getOptionValues("vcf"));
		for (int i = 0; i < vcffiles.length; i++){
			final File vcffile = vcffiles[i];
			if (!vcffile.canRead() || ! vcffile.isFile()){
				log.log(Level.SEVERE, "Can't read vcf file: " + vcffile);
				System.exit(10);
			}
		}
		File[]		bamfiles = Utils.filesFromStrings(line.getOptionValues("bam"));
		String[] 	names;
		if (line.hasOption("samplename")){
			names = line.getOptionValues("samplename");
		} else {
			names = new String[bamfiles.length];
			for (int i = 0; i < bamfiles.length; i++){
				names[i] = bamfiles[i].getName().replace(".bam", "");
			}
		}
		if (names.length != bamfiles.length){
			throw new IllegalArgumentException("Names and BAM files don't match");
		}
		
		File fastafile = new File(line.getOptionValue("r"));
		//IndexedFastaSequenceFile fastaref = new IndexedFastaSequenceFile(new File(line.getOptionValue("r")));

		List<BAMInterface> baminterfaces = new ArrayList<BAMInterface>(names.length);
		for (int i = 0; i < names.length; i++){
			final File bamfile = bamfiles[i];
			if (!bamfile.canRead() || !bamfile.isFile()){
				log.log(Level.SEVERE, "Can't read bam file: " + bamfile);
				System.exit(10);
			}
			final String name = names[i];
			log.log(Level.ALL, "Pairing " + name + " with " + bamfile);
			baminterfaces.add(new BAMInterface(bamfile, name, null));
		}
		// convert to files

		Integer		buffer = Integer.decode(line.getOptionValue("b", "20"));
		Integer		padding = Integer.decode(line.getOptionValue("p", "0"));
		Settings.threadCount = Integer.decode(line.getOptionValue("nThreads", "6"));
		String		outputVCF = line.getOptionValue("o", "wheeljack-alleleconsolidated.vcf");
		Integer		f			=	Integer.decode(line.getOptionValue("f", "0"));
		Integer		F			=	Integer.decode(line.getOptionValue("F", "1284"));
		AlleleResolver.minAlleleCount = Integer.decode(line.getOptionValue("minAlleleCount", "1"));
		Integer		maxSize		=	Integer.decode(line.getOptionValue("maxSize", "20"));

		VCFMergeAndAnnotate merger = new VCFMergeAndAnnotate();
		File vcfoutputFile = new File(outputVCF);
		AlleleResolver.ResolutionType resolution = null;
		if (line.hasOption("p")){
			resolution = AlleleResolver.ResolutionType.NOEXPANDING;
		} else if (line.hasOption("minExp")){
			resolution = AlleleResolver.ResolutionType.MINIMALEXPANDING;
		} else {
			resolution = AlleleResolver.ResolutionType.EXPANDING;
		}
		merger.run(baminterfaces, Arrays.asList(vcffiles), buffer, fastafile, vcfoutputFile, resolution, padding, null, f, F, maxSize);
	}

	public void run(List<BAMInterface> baminterfaces, List<File> variantFiles,
			int buffer, File fastafile, File outputFile,
			AlleleResolver.ResolutionType resolution, int padding,
			File sampleInfo, int f, int F, int maxSize) throws Exception {
		List<String> samples = new ArrayList<String>();
		for (BAMInterface bi : baminterfaces) {
			samples.add(bi.getSampleName());
			log.log(Level.FINE, "Adding sample: " + bi.getSampleName());
		}

		StateMonitor monitor = new StateMonitor();
		MultiVCFReader reader = new MultiVCFReader(variantFiles, buffer, maxSize,
				fastafile);
		// VCFWriter writer = new VCFWriter();
		CARNACSampleGenotyper carnacGenotyper = new CARNACSampleGenotyper();
		VCFHeader vcfHeader = new VCFHeader(
				CARNACSampleGenotyper.getHeaderLines(samples, sampleInfo),
				samples);
		IndexedFastaSequenceFile fastaref = new IndexedFastaSequenceFile(
				fastafile);
		VariantContextWriter writer = new VariantContextWriterBuilder()
				.setReferenceDictionary(fastaref.getSequenceDictionary())
				.setOutputFile(outputFile).build();
		fastaref.close();
		writer.writeHeader(vcfHeader);
		log.log(Level.INFO, "Starting threads");

		int nWorkers = Settings.threadCount > 3 ? Settings.threadCount * 2 : 3; 

		// BlockingQueue<AlleleContainer> container = new LinkedBlockingQueue<AlleleContainer>(
		// 		nWorkers * 2);
		BlockingQueue<AlleleContainer> writercontainer = new LinkedBlockingQueue<AlleleContainer>(
				nWorkers * 2);
		Map<AlleleContainer, VariantContext> resultMap = Collections
				.synchronizedMap(new HashMap<AlleleContainer, VariantContext>());
		// ThreadedAlleleContainerLoader loader = new ThreadedAlleleContainerLoader(
		//		reader, container, writercontainer, monitor);
		// Thread loaderThread = new Thread(loader);
		// loaderThread.start();
		
		// set up the writer
		ThreadedAlleleWriter alleleWriter = new ThreadedAlleleWriter(writer,
				writercontainer, resultMap, monitor);
		Thread writerThread = new Thread(alleleWriter);
		writerThread.start();
		
		
		// using executorservice
		
		final ExecutorService pool = Executors.newFixedThreadPool(nWorkers);
		int allelesAdded = 0;
		while (reader.hasNext()){
			AlleleContainer ac = reader.nextAlleleSet();
			// start executing the thread before the writer is expecting to see it
			pool.execute(new SingleThrowThreadedAlleleResolver(padding, resolution, carnacGenotyper, baminterfaces,
					resultMap, fastafile, ac, monitor, f, F));
			writercontainer.put(ac);
			allelesAdded += 1;
		}
		pool.shutdown();
		while (!pool.isTerminated()){
			try {
				pool.awaitTermination(10, TimeUnit.SECONDS);
				log.fine("Waiting for pool to shut down");
			} catch (InterruptedException e) {
				// this should never happen
				log.fine("Interrupted while shutting down pool.");
			}
		}
		log.log(Level.INFO, "Setting monitor to writing");
		monitor.setWriting();
		// now close the writer thread
		// writerThread.interrupt();
		log.log(Level.INFO, "Joining writer thread");
		writerThread.join();
		writer.close();
		log.log(Level.INFO, "Processed " + allelesAdded + " allele sets.");
		
		/*
		List<Thread> workerThreads = new ArrayList<Thread>(nWorkers);
		workers = new ArrayList<ThreadedAlleleResolver>(nWorkers);
		for (int i = 0; i < nWorkers; i++) {
			final ThreadedAlleleResolver resolver = new ThreadedAlleleResolver(
					padding, resolution, carnacGenotyper, baminterfaces,
					resultMap, fastafile, container, monitor, f, F);
			final Thread resolverThread = new Thread(resolver);
			resolverThread.start();
			workerThreads.add(resolverThread);
			workers.add(resolver);
		}
		
		loaderThread.join();
		log.log(Level.INFO, "Joined the loader thread");
		// set the monitor to processing state.  This means that 
		// the next available thread worker will pull all of the jobs into itself
		// and process them.  This will lead to serial processing, but there won't 
		// be that many sites to work on.
		monitor.setProcessing();
		
		// join into the threads
		for (int i = 0; i < workerThreads.size(); i++){
			final Thread resolverThread = workerThreads.get(i);
			log.log(Level.INFO, "Joining thread " + i + " of " + workerThreads.size());
			resolverThread.interrupt();
			resolverThread.join();
			log.log(Level.INFO, "Thread joined");
		}
		log.log(Level.INFO, "Setting monitor to writing");
		monitor.setWriting();
		// now close the writer thread
		writerThread.interrupt();
		writerThread.join();
		writer.close();
		// finally, collect the new compiled VCF and write to file.
		*/
	}

}
