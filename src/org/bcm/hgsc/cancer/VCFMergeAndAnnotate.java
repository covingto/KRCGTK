package org.bcm.hgsc.cancer;

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
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
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
 * Merge and annotate VCF files.  This merges multiple VCF files together for both tumor and normal samples, extracts relevant sequence information from the BAM files and writes a new VCF file with the information.
 * This method extracts all required information for making a call in the CARNAC variant filter.
 * @author covingto
 *
 */
public class VCFMergeAndAnnotate {
	private static Logger log = Logger.getLogger("");
	private static enum State { READING, PROCESSING, WRITING, DONE};
	
	private class StateMonitor{
		private State state = State.READING;
		StateMonitor(){}
		public void setProcessing(){
			this.state = State.PROCESSING;
		}
		public void setWriting(){
			this.state = State.WRITING;
		}
		public void setDone(){
			this.state = State.DONE;
		}
		public State getState(){
			return this.state;
		}
	}
	private class ThreadedAlleleContainerLoader implements Runnable{
		private final MultiVCFReader reader;
		private final BlockingQueue<AlleleContainer> container;
		private final BlockingQueue<AlleleContainer> writerContainer;
		private boolean done = false;
		private final StateMonitor monitor;

		ThreadedAlleleContainerLoader(MultiVCFReader reader, BlockingQueue<AlleleContainer> container, BlockingQueue<AlleleContainer> writerContainer, StateMonitor monitor){
			this.reader = reader;
			this.container = container;
			this.writerContainer = writerContainer;
			this.monitor = monitor;
		}

		@Override
		public void run(){
			while (this.reader.hasNext() ){
				try {
					AlleleContainer ac = this.reader.nextAlleleSet();
					synchronized (log) {
						log.log(Level.FINE, "Adding allele container " + ac.toString());
					}
					container.put(ac);
					writerContainer.put(ac);
				} catch (InterruptedException e) {
					if (this.reader.hasNext()){
						log.log(Level.SEVERE, "Reader interrupted while adding allele to pool.");
					}
					break;
				}
			}
			this.done = true;
			synchronized(log){
				log.log(Level.INFO, "Done loading allele containers");
			}
			return;
		}

		public boolean done(){
			return done;
		}
	}

	private class ThreadedAlleleWriter implements Runnable{
		private final VariantContextWriter writer;
		private final BlockingQueue<AlleleContainer> alleleContainer;
		private final Map<AlleleContainer, VariantContext> resultMap;
		private boolean done = false;
		private final StateMonitor monitor;

		public ThreadedAlleleWriter(VariantContextWriter writer, BlockingQueue<AlleleContainer> alleleContainer, Map<AlleleContainer, VariantContext> resultMap, StateMonitor monitor){
			this.writer = writer;
			this.alleleContainer = alleleContainer;
			this.resultMap = resultMap;
			this.monitor = monitor;
		}

		@Override
		public void run(){
			while (!Thread.currentThread().isInterrupted() && !monitor.getState().equals(State.WRITING)){
				try {	
					final AlleleContainer ac = this.alleleContainer.take();
					resolve(ac);
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
			}
			final List<AlleleContainer> finalContainers = new ArrayList<AlleleContainer>();
			this.alleleContainer.drainTo(finalContainers);
			for (final AlleleContainer ac : finalContainers){
				resolve(ac);
			} 
			this.done = true;
		}

		private void resolve(AlleleContainer ac){
			while (!resultMap.containsKey(ac)){
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
			}
			final VariantContext vc = resultMap.get(ac);
			if (vc != null){
				writer.add(vc);
			}
			synchronized(resultMap){
				resultMap.remove(ac);
			}
		}
		
		public boolean done(){
			return done;
		}
	}

	private class ThreadedAlleleResolver implements Runnable{
		private final List<BAMInterface> baminterfaces;
		// private final BlockingQueue<VariantContext> variantContextQueue; // will be used only once to add the variant
		// private final AlleleContainer alleleContainer;
		private final AlleleResolver.ResolutionType resolution;
		private final SampleGenotyper genotyper;
		private final File fastafile;
		private final int padding;
		//private final VariantContextWriter writer;
		private final BlockingQueue<AlleleContainer> container;
		private final Map <AlleleContainer, VariantContext> resultMap;
		private final StateMonitor monitor;

		/**
		 * Class initialization should contain all information required to run the method.  Running is just a pass through
		 * @param monitor 
		 */
		public ThreadedAlleleResolver(int padding, AlleleResolver.ResolutionType resolution, 
				SampleGenotyper genotyper, List<BAMInterface> baminterfaces, 
				Map <AlleleContainer, VariantContext> resultMap, File fastafile, BlockingQueue<AlleleContainer> container, StateMonitor monitor){
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
		}

		@Override
		public void run() {
			log.log(Level.FINE, "Starting allele processing");
			while (!Thread.currentThread().isInterrupted() && !monitor.getState().equals(State.PROCESSING)){
				final AlleleContainer alleleContainer = this.container.poll();
				if (alleleContainer != null){
					final VariantContext newv = processAlleleContainer(alleleContainer);
					synchronized(this.resultMap){
						this.resultMap.put(alleleContainer, newv);
					}
				}
			}
			// drain the pool and resolve
			final List<AlleleContainer> finalContainers = new ArrayList<AlleleContainer>();
			synchronized(this.container){
				this.container.drainTo(finalContainers);
			}
			for (final AlleleContainer alleleContainer : finalContainers){
				final VariantContext newv = processAlleleContainer(alleleContainer);
				synchronized(this.resultMap){
					this.resultMap.put(alleleContainer, newv);
				}
			}
		}

		private VariantContext processAlleleContainer(AlleleContainer alleleContainer){
			VariantContext newv = null;
			try{
				final List<ConformedRead> allReads = new ArrayList<ConformedRead>(1000);
				final Map<String, List<ConformedRead>> sampleReads = new HashMap<String, List<ConformedRead>>();
				AlleleSet alleles = null;
				try {
					for (BAMInterface bi : baminterfaces){
						List<ConformedRead> reads = BAMUtils.getConformedReads(
								bi,
								alleleContainer.getChr(), 
								alleleContainer.getStart() - padding, 
								alleleContainer.getEnd() + padding, 
								true,
								new IndexedFastaSequenceFile(this.fastafile)
								);
						//synchronized (log){
						//	log.log(Level.FINE, "Updating reads");
						//}
						sampleReads.put(bi.getSampleName(), reads);
						allReads.addAll(reads);
					}
					//synchronized (log) {
					//	log.log(Level.FINE, "Simplifying allele");
					//	}
					alleles = AlleleResolver.resolveAlleles(allReads, 
							alleleContainer.getChr(), 
							alleleContainer.getStart() - padding, 
							alleleContainer.getEnd() + padding, 
							resolution, new IndexedFastaSequenceFile(this.fastafile), this.padding).simplify();

					// add the info to the variant context builder
					if (alleles == null || alleles.getAlleles().size() == 0){
						// in case the alleles are null, we still want to output something, so we will default to the standard set
						synchronized (log) {
							log.log(Level.FINE, "Adding raw contexts with no genotypes since no alleles returned from allele resolver.");
						}
						for (VariantContext vc : alleleContainer.getVariantContexts()){
							VariantContextBuilder vcbuilder = new VariantContextBuilder();
							vcbuilder.alleles(vc.getAlleles());
							synchronized (log){
								log.log(Level.FINE, "Adding null context to writer");
							}
							newv = vcbuilder.make();
						}

					} else {
						// build the attribute map
						Map<String, Object> varAttributes = new TreeMap<String, Object>();
						varAttributes.put("OC", alleleContainer.getCallString());
						VariantContextBuilder vcbuilder = new VariantContextBuilder();
						
						vcbuilder.alleles(alleles.getAlleles());
						List<Genotype> genotypes = new ArrayList<Genotype>(sampleReads.size());
						// do the genotyping
						for (String k : sampleReads.keySet()){
							final List<ConformedRead> reads = sampleReads.get(k);
							genotypes.add(genotyper.genotype(k, alleles, reads));
						}
						// TODO: filter the genotypes for quality
						vcbuilder.loc(alleleContainer.getChr(), alleles.getStart(), alleles.getEnd())
						.attributes(varAttributes)
						.genotypes(genotypes);
						newv = vcbuilder.make();
					}
				} catch (Exception e){
					log.log(Level.SEVERE, "Hit error in determining genotype: " + alleleContainer.toString(), e);
				} finally {
					// explicitly set items to null this might help the garbage collector, anyway these things are in the young generation anyway...
					allReads.clear();
					for (List<ConformedRead> v : sampleReads.values()){
						v.clear();
					}
					sampleReads.clear();
				}
				log.log(Level.FINE, "Variant resolution complete");
			} catch (Exception e){
				log.log(Level.SEVERE, "Error in resolving allele, will add null to the allele map");
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
		options.addOption("o", true, "output vcf file.  this file will be a vcf formatted file, please add .vcf at the end of the file name [<tumorName>_<normalName>.vcf]");
		options.addOption("h", false, "print help");
		options.addOption("d", false, "Should the program be run in debug mode (outputs all)");
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
		merger.run(baminterfaces, Arrays.asList(vcffiles), buffer, fastafile, vcfoutputFile, resolution, padding, null);
	}

	public void run(List<BAMInterface> baminterfaces, List<File> variantFiles, 
			int buffer, File fastafile, File outputFile, 
			AlleleResolver.ResolutionType resolution, int padding, File sampleInfo) throws Exception{
		List<String> samples = new ArrayList<String>();
		for (BAMInterface bi : baminterfaces){
			samples.add(bi.getSampleName());
			log.log(Level.FINE, "Adding sample: " + bi.getSampleName());
		}
		
		StateMonitor monitor = new StateMonitor();
		MultiVCFReader reader = new MultiVCFReader(variantFiles, buffer, fastafile);
		//VCFWriter writer = new VCFWriter();
		CARNACSampleGenotyper carnacGenotyper = new CARNACSampleGenotyper();
		VCFHeader vcfHeader = new VCFHeader(CARNACSampleGenotyper.getHeaderLines(samples, sampleInfo), samples);
		IndexedFastaSequenceFile fastaref = new IndexedFastaSequenceFile(fastafile);
		VariantContextWriter writer = new VariantContextWriterBuilder().setReferenceDictionary(fastaref.getSequenceDictionary()).setOutputFile(outputFile).build();
		fastaref.close();
		writer.writeHeader(vcfHeader);
		log.log(Level.INFO, "Starting threads");

		int nWorkers = Settings.threadCount > 3 ? Settings.threadCount * 2 : 3; // let's say that there is a super duper amount of IO
		
		// KRC: -2 for the containers is an attempt to throttle the workers, they get to work on things, but can't actually exceed the reader and writer threads capacity without pausing
		BlockingQueue<AlleleContainer> container = new LinkedBlockingQueue<AlleleContainer>(nWorkers - 2);
		BlockingQueue<AlleleContainer> writercontainer = new LinkedBlockingQueue<AlleleContainer>(nWorkers - 2);
		Map <AlleleContainer, VariantContext> resultMap = Collections.synchronizedMap(new HashMap<AlleleContainer, VariantContext>());
		ThreadedAlleleContainerLoader loader = new ThreadedAlleleContainerLoader(reader, container, writercontainer, monitor);
		Thread loaderThread = new Thread(loader);
		loaderThread.start();
		List<Thread> workerThreads = new ArrayList<Thread>(nWorkers);
		for (int i = 0; i < nWorkers; i++){
			final ThreadedAlleleResolver resolver = new ThreadedAlleleResolver(
					padding, resolution, carnacGenotyper, baminterfaces, resultMap, fastafile, container, monitor);
			final Thread resolverThread = new Thread(resolver);
			resolverThread.start();
			workerThreads.add(resolverThread);
		}
		ThreadedAlleleWriter alleleWriter = new ThreadedAlleleWriter(writer, writercontainer, resultMap, monitor);
		Thread writerThread = new Thread(alleleWriter);
		writerThread.start();
		loaderThread.join();
		// now wait for the container (that contains alleles from the writer to drain
		while (container.size() > 0){
			Thread.sleep(1000);
		}
		monitor.setProcessing();
		for (final Thread resolverThread : workerThreads){
			resolverThread.join();
		}
		monitor.setWriting();
		// now close the writer thread
		writerThread.join();
		writer.close();
		// finally, collect the new compiled VCF and write to file.

	}

}
