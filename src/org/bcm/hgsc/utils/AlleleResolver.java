package org.bcm.hgsc.utils;


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.utils.BAMUtils.ConformedRead;


/**
 * Class for resolving alleles at specific loci given a set of bam files and input coordinates.
 * Initially candidate regions can be obtained from a VCF or more generally a VariantContext object.
 * However, in some cases alleles may be missed or ambiguous within the context because of alignment or caller differences.
 * 
 * This class will attempt to identify all reasonable alleles within a given loci.
 * @author covingto
 *
 */
public class AlleleResolver {
	private static Logger log = Logger.getLogger(AlleleResolver.class.getName());
	public static enum ResolutionType { 
		EXPANDING, MINIMALEXPANDING, NOEXPANDING
	}
	public static int minAlleleCount = 2;
	public static class AlleleSet{
		private static Logger log = Logger.getLogger(AlleleSet.class.getName());
		private final Set<Allele> alleles;
		private final byte[] reference;
		private final int start;
		private final int end;
		private final int sliceStart;
		private final int sliceEnd;
		private final String contig;
		private final int buffer;

		public AlleleSet(String contig, int start, int end, Collection<Allele> alleles, byte[] reference, int buffer) throws Exception{
			this(contig, start, end, alleles, reference, start, end, buffer);
		}

		public AlleleSet(String contig, int start, int end, Collection<Allele> alleles, byte[] reference, int sliceStart, int sliceEnd, int buffer) throws Exception{
			this.start = start;
			this.sliceStart = sliceStart;
			this.end = end;
			this.sliceEnd = sliceEnd;
			this.alleles = new HashSet<Allele>(alleles);
			this.reference = reference;
			this.contig = contig;
			this.buffer = buffer;
			if (this.reference == null){
				throw new Exception ("reference may not be null");
			} else {
				// add the reference to the set anyway, note that the alleles may already contain the reference, that's OK.
				this.alleles.add(Allele.create(this.reference, true));
			}
		}
		
		@Override
		public String toString(){
			return this.contig + ":" + this.start + "(" + this.sliceStart + ")" + "-" + this.end + "(" + this.sliceEnd + ")" + "\n" +
					StringUtils.join(this.alleles, ", ") + "[" + new String(this.reference) + "]";
		}

		public int getSliceStart(){
			return this.sliceStart;
		}

		public int getSliceEnd(){
			return this.sliceEnd;
		}

		public int getStart(){
			return this.start;
		}

		public int getEnd(){
			return this.end;
		}


		public Collection<Allele> getAlleles() {
			// returns an effective copy of the allele list.  Hopefully this won't come back to bite us. :)
			return this.alleles;
		}

		/**
		 * simplifies and allele set and returns a minimal set of alleles.  All positions of alleles are shifted to the right until
		 * either i-1 position if i is an indel or i if i is a snp
		 * 
		 * From the example in resolveAlleles;
		 * Alleles all beginning at position 6 (the minimal position available) in the form;
		 * AGGGTTCCC	(reference)
		 * AGGCTCCC		(read2)
		 * AGGCTCCC		(read3, note that this is the SAME allele as in read 2)
		 * AGGCTCAC		(read5)
		 * AGGGTTCCCC	(read6)
		 * 
		 * returns;
		 * position 8
		 * GGTTCCC
		 * GCTCCC
		 * GCTCAC
		 * GGTTCCCC
		 * @return
		 */
		public AlleleSet simplify() {
			try {
				// check for the simple and most common case where there is only a SNP
				boolean isSNP = false; // initialized to false, we will set this to true if the alleles only support SNP
				int maxLen = reference.length;
				int minLen = Integer.MAX_VALUE;
				int rightoffset = Integer.MAX_VALUE;
				int leftoffset = Integer.MAX_VALUE;
				// fallback if there are no alleles
				if (this.alleles.size() == 0 ){ 
					log.finer("Allele length is 0");
					return this;
				}
				for (Allele a : this.alleles){
					//log.log(Level.FINE, "Checking " + a.toString());
					maxLen = maxLen > a.length() ? maxLen : a.length();
					minLen = minLen < a.length() ? minLen : a.length();
					if (!a.isReference()){
						
						final int a_rightoffset = getRightOffset(a.getBases(), reference);
						final int a_leftoffset = correctLeftOffset(a.getBases(), reference, getLeftOffset(a.getBases(), reference), a_rightoffset);
						// get the min of each left and right offset
						rightoffset = a_rightoffset < rightoffset ? a_rightoffset : rightoffset;
						leftoffset = a_leftoffset < leftoffset ? a_leftoffset : leftoffset;
					}
				}
				if ((maxLen == 1) || (leftoffset == 0 && rightoffset == 0)){
					log.log(Level.FINER, "No reduction required");
					return this; // this is the simplest way in which to represent these data
				}
				
				// special case: if minLen == maxLen then this is actually an onp (snp, dnp, tnp, etc)
				// therefore we must index up the leftoffset 
				if (minLen == maxLen){
					isSNP = true;
				} else {
					isSNP = false;
					leftoffset = 0 > (leftoffset - 1) ? 0 : (leftoffset - 1); // decriment the left offset since this is an indel it must be anchored by reference
				}
				
				if (leftoffset >= (minLen - rightoffset)){
					log.log(Level.FINER, "Reduction not possible: minLength: " + minLen + " leftoffset: " + leftoffset + " rightoffset: " + rightoffset + "\nAlleles: " + StringUtils.join(this.alleles, ", "));
					return this;
				}
				
				//
				// This section uses the offsets calculated from above and generates the new allele sets that will be used, removing those that are the same as reference.
				
				Set<Allele> newAlleles = new HashSet<Allele>();
				// Filter alleles from the set, we are just awash in alleles that are really of low quality.  Added a freature that an allele must be present at least a minimal number of times to be considered valid.
				try {
					// handle the reference allele
					byte[] newreference = this.reference;
					// if (leftoffset < this.reference.length - rightoffset){
					newreference = Arrays.copyOfRange(this.reference, leftoffset, this.reference.length - rightoffset);
					//} else {
						// this is very strange, what alleles were sequenced anyway?
						// no matter, clear the offsets and just go with what we have
					//	log.warning("Left and right offsets do not conform: " + leftoffset + " !< " + this.reference.length + " - " + rightoffset);
					//	leftoffset = 0;
					//	rightoffset = 0;
					//}
					// the reference is always added
					final Allele newReferenceAllele = Allele.create(newreference, true);
					newAlleles.add(newReferenceAllele);
					log.finest("Reference bases are: " + newReferenceAllele.getBaseString());
					for (Allele a : this.alleles){
						final byte[] tmpAllele = Arrays.copyOfRange(a.getBases(), leftoffset, a.length() - rightoffset);
						try {
							final Allele newAllele = Allele.create(tmpAllele, false);
							// only add if the allele does not match the reference, others are handled in the set by .equals in the set
							if (!newReferenceAllele.basesMatch(newAllele)){ 
								log.finest("Plan to add allele with bases: " + newAllele.getBaseString());
								// this is a non-reference allele, so we get to add it if it passes some checks
								if (isSNP){
									if (tmpAllele.length == newreference.length){
										newAlleles.add(newAllele);
									} else {
										log.warning("Wanted a SNP but allele length was not 1, this really shouldn't happen");
									}
								} else {
									newAlleles.add(newAllele);
								}
							}
						} catch (Exception e) {
							log.log(Level.SEVERE, "Error in comparing alleles: " + a.getBaseString() + "; reference = " + new String(this.reference) + " leftoffset = " + leftoffset + " len - rightoffset = " + (a.length() - rightoffset));
							throw e;
						}
					}
					
				} catch (Exception e) {
					//System.out.println("Error after offset calculations; leftoffset: " + leftoffset + " rightoffset: " + rightoffset);
					log.log(Level.SEVERE, e.toString(), e);
					throw e;
				}
				log.finer("leftoffset: " + leftoffset + " rightoffset: " + rightoffset);
				return new AlleleSet(this.contig, this.start + leftoffset, this.end - rightoffset, newAlleles, Arrays.copyOfRange(this.reference, leftoffset, this.reference.length - rightoffset), this.start, this.end, buffer);
			} catch (Exception e){
				log.log(Level.WARNING, e.toString(), e);
				//e.printStackTrace();
				log.log(Level.WARNING, "Printing alleles relating to error");
				for (Allele a : this.alleles){
					log.log(Level.WARNING, a.toString());
				}
				return this;
			}
		}
		
		public static int getRightOffset(byte[] a, byte[] b){
			int minLength = a.length < b.length ? a.length : b.length;
			for (int i = 0; i < minLength; i++){
				if (getByteAtRightOffset(a, i) != getByteAtRightOffset(b, i)){
					return i;
				}
			}
			return 0; // default is no offsetting, this is less risky
		}
		
		public static int getLeftOffset(byte[] a, byte[] b){
			int minLength = a.length < b.length ? a.length : b.length;
			for (int i = 0; i < minLength; i++){
				if (getByteAtLeftOffset(a, i) != getByteAtLeftOffset(b, i)){
					return i;
				}
			}
			return minLength - 1;
		}
		
		public static int correctLeftOffset(byte[] a, byte[] b, int left, int right){
			int minlen = a.length < b.length ? a.length : b.length;
			if ((left + right) > minlen){
				int newleft = left - (minlen - left - right) * -1;
				return newleft;
			} else {
				return left;
			}
		}
		
		// sure this could be in-lined but why not keep it clear
		private static byte getByteAtLeftOffset(byte[] ba, int offset) {
			// offset 0 ==> return first base
			// offset ba.length ==> return last base
			if (offset >= ba.length){
				throw new IndexOutOfBoundsException("Offset is greater than the array length");
			}
			return ba[offset];
		}

		private static byte getByteAtRightOffset(byte[] ba, int offset){
			// offset 0 ==> return last base
			// offset ba.length ==> return first base
			if (offset >= ba.length){
				throw new IndexOutOfBoundsException("Offset is greater than the array length");
			}
			return ba[ba.length - (offset + 1)];
		}

		public String getChr() {
			return this.contig;
		}
	}

	/**
	 * Pass through the available reads and generate the start and end locations of the alleles.  
	 * These may later be paired down and simplified in later steps (AlleleSet.simplify()).
	 * 
	 * Pos				[	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14		15	16	17	]
	 * Reference:		[	A	T	C	T	C	T	A	G	G	G	T	T	C	C	C	-	T	G	G	]
	 * Read1:			[	A	T	C	T	C	T	A	G	G	G	T	T	C	C	C	-	T	G	G	]	Pure reference
	 * Read2:			[	A	T	C	T	C	T	A	-	G	G	C	T	C	C	C	-	T	G	G	]	6; AGGGT > AGGC
	 * Read3:			[	A	T	C	T	C	T	A	G	G	C	-	T	C	C	C	-	T	G	G	]   8; GGT > GC
	 * Read4:			[	A	T	C	T	C	T	A	G	G	G	T	T	C	C	C	-	T	G	G	]	Pure reference
	 * Read5:			[	A	T	C	T	C	T	A	-	G	G	C	T	C	A	C	-	T	G	G	]	6; AGGGTTCC > AGGCTCA
	 * Read6:			[	A	T	C	T	C	T	A	G	G	G	T	T	C	C	C	C	T	G	G	]	14; C > CC
	 * 
	 * These will result in the generation of bounded alleles all beginning at position 6 (the minimal position available) in the form;
	 * AGGGTTCCC	(reference)
	 * AGGCTCCC		(read2)
	 * AGGCTCCC		(read3, note that this is the SAME allele as in read 2)
	 * AGGCTCAC		(read5)
	 * AGGGTTCCCC	(read6)
	 * 
	 * @param reads
	 * @param start
	 * @param end
	 * @param ignoreExtra
	 * @return
	 * @throws Exception 
	 */
	public static AlleleSet resolveAlleles(List<ConformedRead> reads, String contig, int start, int end, ResolutionType resolution, IndexedFastaSequenceFile fastaref, int buffer) throws Exception{
		log.log(Level.FINEST, "Processing " + reads.size() + " reads");
		List<Allele> alleles = new ArrayList<Allele>();
		int cri = 0;
		parseReads: for (final ConformedRead cr : reads){
			cri += 1;
			if (cr.readStart() > start || cr.readEnd() < end){ continue; }
			//System.out.println("Debug: parsing conformed read");
			//System.out.println(cr.toString());
			final int[] crPos = cr.getAlleleRangeAtGenomicPos(start);
			final int crStart = crPos[0];
			final int crEnd = crPos[1];
			//if (crEnd - crStart > 50){ continue; } // don't expand alleles by 50, these will be independend records

			// this is the case where we need to expand the allele positions.  New positions will be included in the new AlleleSet returned.
			switch (resolution){
				case EXPANDING:
					if (crStart < start && crEnd - crStart < 50){
						return resolveAlleles(reads, contig, crStart, end, resolution, fastaref, buffer);
					} else if (crEnd > end && crEnd - crStart < 50){
						return resolveAlleles(reads, contig, start, crEnd, resolution, fastaref, buffer);
					} else {
						// There is no need to expand here
						break;
					}
				case MINIMALEXPANDING:
					if (crStart < start){
						return resolveAlleles(reads, contig, crStart, end, resolution, fastaref, buffer);
					} else {
						// There is no need to expand here
						break;
					}
				case NOEXPANDING:
					break;
//					if (crStart < start || crEnd > end){
//						log.log(Level.FINEST, "Missed expansion of allele " + cr.toString() + " crStart: " + crStart + " start: " + start + " crEnd: " + crEnd + " end: " + end);
//						//continue parseReads;
//					} else {
//						// There is no need to expand here
//						//log.log(Level.WARNING, "Missed expansion of allele " + cr.toString() + " crStart: " + crStart + " start: " + start + " crEnd: " + crEnd + " end: " + end + " breaking!!!");
//						break;
//					}
				default:
					{
						log.log( Level.WARNING, "Hit theoretically unreachable code");
						continue parseReads;
					}
			}
			// create a new Allele by slicing the conformed read to start and end
			byte[] seqAllele = cr.getSeqAllele(start, end).bytes;
			// log.log(Level.FINEST, "Checking allele " + new String(seqAllele));
			if (seqAllele.length < 1){
				log.log(Level.SEVERE, "Generated empty allele" + cr.toString() + " crStart: " + crStart + " start: " + start + " crEnd: " + crEnd + " end: " + end);
				continue parseReads;
			}
			for (int i = 0; i < seqAllele.length; i++){
				if (seqAllele[i] == BAMUtils.dot || seqAllele[i] == BAMUtils.n || seqAllele[i] == BAMUtils.N|| seqAllele[i] == BAMUtils.zero || seqAllele[i] == BAMUtils.unk){
					log.log(Level.FINEST, "Nonconforming base at position i = " + i + " in " + new String(seqAllele));
					continue parseReads;
				}
			}
			// log.log(Level.FINEST, "Processed " + cri);
			alleles.add(Allele.create(seqAllele, false));
		}
		log.log(Level.FINE, "Generated allele set with " + alleles.size() + " acceptable reads");
		
		// initial simplification of the allele set
		// now that the variants are added to the condensedAlleles, we count the alleles.
		Set<Allele> newAlleles = new HashSet<Allele>();
		Map<Allele, Integer> alleleCount = Utils.countOccurrences(alleles);
		// filter the alleles
		final Allele referenceAllele = Allele.create(SynchronousIndexedFastaReader.getSubsequenceAt(fastaref, contig, start, end).getBases(), true);
		newAlleles.add(referenceAllele);
		for (final Entry<Allele, Integer> entry : alleleCount.entrySet()){
			if (referenceAllele.basesMatch(entry.getKey())){
				continue;
			}
			if (entry.getValue() >= AlleleResolver.minAlleleCount){
				newAlleles.add(entry.getKey());
			} 
			else {
				log.log(Level.FINEST, "Discarded allele " + entry.getKey().toString() + " because of insufficient coverage (" + entry.getValue() + "<" + AlleleResolver.minAlleleCount + "). Actual frequency was " + Collections.frequency(alleles, entry.getKey()));
			}
		}
		return new AlleleSet(contig, start, end, newAlleles, SynchronousIndexedFastaReader.getSubsequenceAt(fastaref, contig, start, end).getBases(), buffer);
	}
}
