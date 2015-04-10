package org.bcm.hgsc.utils;


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

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
				// add the reference to the set anyway
				this.alleles.add(Allele.create(this.reference, true));
			}
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
				// fallback if there are no alleles
				if (this.alleles.size() == 0 ){ 
					return this;
				}
				for (Allele a : this.alleles){
					//log.log(Level.FINE, "Checking " + a.toString());
					maxLen = maxLen > a.length() ? maxLen : a.length();
					minLen = minLen < a.length() ? minLen : a.length();
				}
				if (maxLen == 1){
					//log.log(Level.FINER, "No reduction required");
					return this; // this is the simplest way in which to represent these data
				}
				//log.log(Level.FINE, "minlength: " + minLen + " maxlength: " + maxLen);
				// pass through all alleles in the set to get the right chewback amount
				int rightoffset = 0;
				findrightposition: for (int i = 0; i < minLen; i++){
					for (Allele a : this.alleles){
						if (a.isReference()){
							continue;
						}
						//System.out.println(a.toString() + " length: " + a.length() + " chewback: " + i);
						//System.out.println(new String(reference) + " length: " + reference.length + " chewback: " + i);
						// add a -1 offset to the length
						if (a.getBases()[a.length() - i - 1] != reference[reference.length - i - 1]){
							break findrightposition;
						}
					}
					rightoffset = i + 1;
				}
				
				// pass through all alleles in the set to get the last (rightmost) agreeing position
				// we must have checked the rightoffset first as this will define where the alleles must end
				int leftoffset = 0;
				// NOTE: this.buffer > 0 indicates that the alleles were initially buffered, so we can't collapse beyond the base size.
				findleftposition: for (int i = (this.buffer > 0 && minLen > 1 ? 1 : 0); i < minLen - rightoffset; i++){
					for (Allele a : this.alleles){
						if (a.isReference()){
							continue;
						}
						if (a.getBases()[i] != reference[i]){
							break findleftposition;
						}
					}
					leftoffset = i; // this is never reached if the second base is not all equal
				}
				
				// special case: if minLen == maxLen then this is actually an onp (snp, dnp, tnp, etc)
				// therefore we must index up the leftoffset 
				if (minLen == maxLen){
					leftoffset++;
					isSNP = true;
				}
				
				Set<Allele> newAlleles = new HashSet<Allele>();
				// Filter alleles from the set, we are just awash in alleles that are really of low quality.  Added a freature that an allele must be present at least a minimal number of times to be considered valid.
				try {
					// handle the reference allele
					byte[] newreference = this.reference;
					if (leftoffset < this.reference.length - rightoffset){
						newreference = Arrays.copyOfRange(this.reference, leftoffset, this.reference.length - rightoffset);
					} else {
						// this is very strange, what alleles were sequenced anyway?
						// no matter, clear the offsets and just go with what we have
						log.warning("Left and right offsets do not conform: " + leftoffset + " !< " + this.reference.length + " - " + rightoffset);
						leftoffset = 0;
						rightoffset = 0;
					}
					// the reference is always added
					newAlleles.add(Allele.create(newreference, true));
					for (Allele a : this.alleles){
						byte[] tmpAllele = Arrays.copyOfRange(a.getBases(), leftoffset, a.length() - rightoffset);
						// tmpAllele[0] == newreference[0] ensures that this isn't a runaway alignment.  The reference and position are anchored by a similar base
						try {
							if (! Arrays.equals(tmpAllele, newreference)){ 
								// this is a non-reference allele, so we get to add it if it passes some checks
								Allele newAllele = Allele.create(tmpAllele, false);
								if (isSNP){
									if (tmpAllele.length == newreference.length){
										
										newAlleles.add(newAllele);
									} else {
										log.warning("Wanted a SNP but allele length was not 1, this really shouldn't happen");
									}
								} else {
									// Commented by KRC; issue is when the variant actually is an indel but there are variant reads to the left and right
									//if (tmpAllele[0] != newreference[0]){
									//	log.warning("Wanted an indel but anchoring base is not correct: " + a.toString());
									//} else {
									newAlleles.add(newAllele);
									//}
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

		List<Allele> alleles = new ArrayList<Allele>();
		parseReads: for (ConformedRead cr : reads){
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
						break;
					}
				case MINIMALEXPANDING:
					if (crStart < start){
						return resolveAlleles(reads, contig, crStart, end, resolution, fastaref, buffer);
					} else {
						break;
					}
				case NOEXPANDING:
					if (crStart < start || crEnd > end){
						// log.log(Level.WARNING, "Missed expansion of allele " + cr.toString() + " crStart: " + crStart + " start: " + start + " crEnd: " + crEnd + " end: " + end);
						continue parseReads;
					} else {
						break;
					}
				default:
					log.log( Level.WARNING, "Hit theoretically unreachable code");
					continue parseReads;
			}
			// create a new Allele by slicing the conformed read to start and end
			byte[] seqAllele = cr.getSeqAllele(start, end).bytes;
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
//			checkReference: if (reference == null){
//				byte[] tmpreference = ConformedRead.removePlaceholders(cr.getRefAtGenomicRange(start, end)).bytes;
//				for (int i = 0; i < tmpreference.length; i++){
//					if (tmpreference[i] == BAMUtils.dot || tmpreference[i] == BAMUtils.n || tmpreference[i] == BAMUtils.zero || tmpreference[i] == BAMUtils.unk){
//						log.log(Level.ALL, "Reference lookup brokd at i = " + i + " tmpreference = " + new String(tmpreference));
//						break checkReference;
//					}
//				}
//				reference = tmpreference;
//				alleles.add(Allele.create(reference, true));
//			}
			
			// the reference sequence is questionable around soft clipped bases, therefore, we must compare the allele at this stage and only add if the allele is explicitly reference
			//System.out.println("Adding allele to set");
			alleles.add(Allele.create(seqAllele, false));
		}
		
		//System.out.println("Allele set complete");
//		if (reference == null){
//			throw new Exception("Reference can not be null");
//		}
		log.log(Level.FINE, "Generated allele set with " + alleles.size() + " acceptable reads");
		
		// initial simplification of the allele set
		// now that the variants are added to the condensedAlleles, we count the alleles.
		Set<Allele> newAlleles = new HashSet<Allele>();
		Map<Allele, Integer> alleleCount = Utils.countOccurrences(alleles);
		// filter the alleles
		for (Entry<Allele, Integer> entry : alleleCount.entrySet()){
			if (entry.getValue() >= AlleleResolver.minAlleleCount){
				newAlleles.add(entry.getKey());
			} 
			//else {
			//	log.log(Level.FINE, "Discarded allele " + entry.getKey().toString() + " because of insufficient coverage " + entry.getValue() + " Actual frequency was " + Collections.frequency(alleles, entry.getKey()));
			//}
		}
		return new AlleleSet(contig, start, end, newAlleles, fastaref.getSubsequenceAt(contig, start, end).getBases(), buffer);
	}
}
