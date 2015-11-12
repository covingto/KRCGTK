package org.bcm.hgsc.cancer;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.lang3.StringUtils;


/**
 * Reads data from multiple VCF files simultaneously to generate allele specific data.
 * Has the ability to return and resolve multiple alleles given a buffer size.
 * 
 * Internally this opens multiple streams to each of the input files for each of the different samples.  These are then read sequentially and returned 
 * in sets based on the buffer size requested.  For a simple case imagine the following;
 * 
 * Sample 1:
 * 		Indel file:
 * 			chr		pos		ref		var
 * 			1		10		AC		A
 * 			1		20		G		GT
 * 			1		500		ATT		A
 * 
 * 		SNP file:
 * 			1		15		A		T
 * 			1		300		T		G
 * 
 * Sample 2:
 * 		Indel file
 * 			1		10		AC		A
 * 			1		500		ATT		A
 * 
 * 		SNP file:
 * 			1		15		A		T
 * 
 * With a buffer size of 50 this will return the following allele containers:
 * 	Allele Container
 * 		variants:		[(1,10,AC,A),(1,10,AC,A),(1,15,A,T),(1,15,A,T),(1,20,G,GT)]
 * 		alleles:		[(1,10,AC),(1,10,A),(1,10,AC),(1,10,A),(1,15,A),(1,15,T),...]
 * 		callerPaths:	[sample1/indelfile, sample2/indelfile, sample1/snpfile, sample2/snpfile, sample1/indelfile]
 * 		start:			10
 * 		end:			20
 * 		chr:			1
 * 
 * 	Allele Container
 * 		variants:		[(1,300,T,G)]
 * 		alleles:		[(1,300,T),(1,300,G)]
 * 		callerPaths:	[sample1/snpfile]
 * 		start:			300
 * 		end:			300
 * 		chr:			1
 * 
 * Allele Container
 * 		variants:		[(1,500,ATT,A),(1,500,ATT,A)]
 * 		alleles:		[(1,500,ATT),(1,500,A),(1,500,ATT),(1,500,A)]
 * 		callerPaths:	[sample1/indelfile, sample2/indelfile]
 * 		start:			500
 * 		end:			500
 * 		chr:			1
 *
 * @author covingto
 *
 */
public class MultiVCFReader {
	private static Logger log = Logger.getLogger(MultiVCFReader.class.getName());
	private final int buffer; // the buffer size for finding variants
	//private final Map<String, CloseableTribbleIterator<VariantContext>> tumorReaders = new HashMap<String, CloseableTribbleIterator<VariantContext>>();
	//private final Map<String, CloseableTribbleIterator<VariantContext>> normalReaders = new HashMap<String, CloseableTribbleIterator<VariantContext>>();
	private final Map<String, CloseableTribbleIterator<VariantContext>> vcfReaders = new HashMap<String, CloseableTribbleIterator<VariantContext>>();
	//private final Map<String, VariantContext> tumorNext = new HashMap<String, VariantContext>();
	//private final Map<String, VariantContext> normalNext = new HashMap<String, VariantContext>();
	private final Map<String, VariantContext> vcfNext = new HashMap<String, VariantContext>();
	private final SAMSequenceDictionary sequenceDict;
	private SAMSequenceRecord currentSequenceRegion = null;
	private int currentSequenceRecordIndex = 0;
	private String lastSource;
	private VariantContext nextVariantContext;
	private final int maxSize;
	
	static {
		log.setLevel(Level.ALL);
	}
	
	public class AlleleContainer{
		private final List<VariantContext> variants;
		private final List<Allele> alleles;
		private final List<String> callerPaths;
		private final int start;
		private final int end;
		private final String chr;
		
		AlleleContainer(String chr, int start, int end, List<VariantContext> variants, List<String> callerPaths, List<Allele> alleles){
			this.variants = variants;
			this.alleles = alleles;
			this.callerPaths = callerPaths;
			this.chr = chr;
			this.start = start;
			this.end = end;
		}
		
		@Override
		public String toString(){
			return "Start: " + this.start + " End: " + this.end + " Alleles: " + this.alleles.size();
		}
		
		public String getChr(){
			return this.chr;
		}
		
		public int getStart(){
			return this.start;
		}
		
		public int getEnd(){
			return this.end;
		}
		
		public List<VariantContext> getVariantContexts(){
			return new ArrayList<VariantContext>(this.variants);
		}
		
		public List<Allele> getAlleles(){
			return new ArrayList<Allele>(this.alleles);
		}
		
		public List<String> getCallerPaths(){
			return new ArrayList<String>(this.callerPaths);
		}
		
		public String getCallString(){
			List<VariantContext> variants = this.getVariantContexts();
			List<String> callerPaths = this.getCallerPaths();
			// check that the variants paths and variants are the same length
			String[] callerHits = new String[callerPaths.size()];
			for (int i = 0; i < variants.size(); i++){
				final VariantContext thisVariant = variants.get(i);
				final String thisCallerPath = callerPaths.get(i);
				final List<Allele> variantAlleles = thisVariant.getAlternateAlleles();
				StringBuilder sbuilder = new StringBuilder(thisCallerPath + "(");
				for (int j = 0; j < variantAlleles.size(); j++){
					final Allele ref = thisVariant.getReference();
					final Allele a = variantAlleles.get(j);
					final int pos = thisVariant.getStart();
					sbuilder.append("[" + pos + "~");
					sbuilder.append(ref.getBaseString() + ">");
					sbuilder.append(a.getBaseString() + "]");
				}
				sbuilder.append(")");
				callerHits[i] = sbuilder.toString();
			}
			return StringUtils.join(callerHits, "|");
		}
		
		
	}
	
	MultiVCFReader(List<File> vcfFiles, int buffer, int maxSize, File fastaFile) throws FileNotFoundException, Exception{
		this(vcfFiles, buffer, maxSize, new IndexedFastaSequenceFile(fastaFile));
	}
	
	MultiVCFReader(List<File> vcfFiles, int buffer, int maxSize, IndexedFastaSequenceFile fastaref) throws Exception{
		this.buffer = buffer;
		this.sequenceDict = fastaref.getSequenceDictionary();
		this.maxSize = maxSize;
		if (this.sequenceDict == null){
			throw new Exception("Fasta sequence dictionary is null");
		}
		this.currentSequenceRegion = this.sequenceDict.getSequence(this.currentSequenceRecordIndex);
		for (File vcf : vcfFiles){
			final VCFCodec vcfCodec = new VCFCodec();
			boolean requireIndex=false;
			FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(
					vcf.getAbsolutePath(), vcfCodec, requireIndex);
			final String key = vcf.getAbsolutePath();
			vcfReaders.put(key, reader.iterator());
			vcfNext.put(key, vcfReaders.get(key).next());
		}
		this._setNewNext(); // set up the iterator
		//try {
			// log.log(Level.ALL, "Next context: " + this.nextVariantContext.toString());
		//} catch (Exception e) {
		//	log.log(Level.ALL, "Caught exception from file: " + this.lastSource);
		//	throw e;
		//}
		
	}
	
	private VariantContext _next(){
		return this.nextVariantContext;
	}
	
	private VariantContext _takeNextFromReader(String source){
		return _takeNextFromReader(vcfReaders.get(source));
	}
	
	private VariantContext _takeNextFromReader(CloseableTribbleIterator<VariantContext> reader){
		try {
			while (reader.hasNext()){
				VariantContext next = reader.next();
				if (next.getEnd() - next.getStart() > this.maxSize){ continue; }
				return next;
			}
			return null;
		} catch (Exception e){
			log.log(Level.SEVERE, "Unacceptable record returned from VCF read from source: " + reader.toString(), e);
			return null;
		}
	}
	
	/** Take the variant in next and index up the variants */
	private VariantContext _takeNext(){
		VariantContext oldNext = this._next();
		if (vcfReaders.containsKey(this.lastSource)){
			if (vcfReaders.get(this.lastSource).hasNext()){
				// TODO: What should happen when next fails because of malformed lines?
				VariantContext c = null;
				final CloseableTribbleIterator<VariantContext> reader = vcfReaders.get(this.lastSource);
				while (reader.hasNext() && c == null){
					c = _takeNextFromReader(this.lastSource); // take may throw errors that are captured, but the reader doesn't stop
				}
				if (c == null){
					vcfNext.remove(this.lastSource);
				} else {
					vcfNext.put(this.lastSource, c);
				}
			} else {
				vcfNext.remove(this.lastSource);
			}
		}
		// index up the next container that was the source for the last next
//		if (vcfNext.containsKey(this.lastSource)){
//			vcfNext.put(this.lastSource, vcfReaders.get(this.lastSource).next());
//		} else {
//			normalNext.put(this.lastSource, normalReaders.get(this.lastSource).next());
//		}
		
		// set the new next
		this._setNewNext();
		
		return oldNext;
	}
	
	private void _setNewNext() {
		VariantContext nvc = null;
		String nsource = null;
		
		for (String k : vcfNext.keySet()){
			VariantContext vc = vcfNext.get(k);
			if (vc == null){ continue; }
			if (!vc.getChr().equals(this.currentSequenceRegion.getSequenceName())){ continue; } // can't add since this is not the current record
			if (nvc == null){
				nvc = vc;
				nsource = k;
			} else {
				if (vc.getStart() < nvc.getStart()){
					nvc = vc;
					nsource = k;
				}
			}
		}

		// now we have the minimal set (we hope)
		// check if the nvc is null, this will happen when we move chromosomes (~24 times in a run) If that is the case, index the chr and update
		if (nvc == null){
			this.currentSequenceRecordIndex++;
			// log.log(Level.FINEST, "Moving to next chromosome contig: " + this.currentSequenceRecordIndex);
			// System.out.println("Moving to next chromosome contig: " + this.currentSequenceRecordIndex);
			// System.out.println("seqdictSize: " + this.sequenceDict.size() + " currentindex: " + this.currentSequenceRecordIndex);
			if (this.currentSequenceRecordIndex < this.sequenceDict.size()){
				this.currentSequenceRegion = this.sequenceDict.getSequence(this.currentSequenceRecordIndex);
				log.log(Level.INFO, "New sequence region: " + this.currentSequenceRegion.getSequenceName());
				// System.out.println("New sequence region: " + this.currentSequenceRegion);
				this._setNewNext();
			} else {
				log.log(Level.INFO, "End of readers");
				// System.out.println("End of readers");
				this.nextVariantContext = null;
				this.lastSource = null;
			}
		} else {
			this.nextVariantContext = nvc;
			this.lastSource = nsource;
		}
	}

	public AlleleContainer nextAlleleSet(){
		log.log(Level.FINE, "Collecting allele set");
		List<Allele> alleleSet = new LinkedList<Allele>(); // add the new alleles to this.
		List<VariantContext> contexts = new LinkedList<VariantContext>(); // the variant contexts from where the alleles will be derived
		List<String> initialCallPaths = new LinkedList<String>(); // collects the caller paths by checking lastSource
		// find the very next allele to grab;
		
		initialCallPaths.add(this.lastSource);
		VariantContext vc = this._takeNext(); // get the next vc in the stack
		if (vc == null){ return null; } // this means we have reached the end, this state will lock here
		contexts.add(vc);
		
		String chr = vc.getChr();
		int contextsStart = vc.getStart();
		//int firstContextEnd = vc.getEnd();
		int contextsEnd = vc.getEnd();
		
		// grab the nearby alleles within the buffer.
		while (this._next() != null && vc.getChr().equals(this._next().getChr()) && contextsEnd > (this._next().getStart() - this.buffer)){
			initialCallPaths.add(this.lastSource);
			vc = this._takeNext();
			contextsEnd = Math.max(contextsEnd, vc.getEnd());
			contexts.add(vc);
			alleleSet.addAll(vc.getAlleles());
		}
		
		return new AlleleContainer(chr, contextsStart, contextsEnd, contexts, initialCallPaths, alleleSet);
		
	}
	
	public boolean hasNext(){
		return this._next() != null;
	}

}
