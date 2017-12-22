package org.bcm.hgsc.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;

/**
 * Utility functions for working with BAM files.
 * @author covingto
 *
 */
public class BAMUtils {
	public static final byte a='a', c='c', g='g', t='t', n='n', A='A', C='C', G='G', T='T', N='N', dot='.', unk='-', zero='0';
	public enum VarType{
		SNV, INS, DEL, ANY
	}
	private static final Logger log = Logger.getLogger(BAMUtils.class.getName());

	public class Pileup{
		private Vector<Byte> seq = new Vector<Byte>();
		private Vector<Byte> qual = new Vector<Byte>();
		private Vector<Integer> orientation = new Vector<Integer>();
		private final byte ref;
		private final String chr;
		private final long pos;
		
		protected Pileup(String chr, long pos, byte ref){
			this.chr = chr;
			this.pos = pos;
			this.ref = ref;
		}
		
		protected void addReadData(byte rs, byte rq, int ro){
			seq.add(rs);
			qual.add(rq);
			orientation.add(ro);
		}
		
		public byte getRef(){
			return this.ref;
		}
		
		public String getChr(){
			return this.chr;
		}
		
		public long getPos(){
			return this.pos;
		}
		
		public Byte[] getSeq(){
			return seq.toArray(new Byte[seq.size()]);
		}
		
		public Byte[] getQual(){
			return qual.toArray(new Byte[qual.size()]);
		}
		
		public Integer[] getOrientation(){
			return orientation.toArray(new Integer[orientation.size()]);
		}
	}

	/**
	 * This class holds a set of arrays representing a conformed read.  This can be used to find variants at positions.
	 * Ex;
	 * 	read;	ATCGATCGATCG 
	 * 	qual;	>>>><<<<>>>>	
	 * 	CIGAR;	3M1I3M1D5M
	 * 	pos;	1
	 * to 
	 * 	ref:	[	A	T	C	-	A	T	C	?	G	A	T	C	G	]
	 * 	read:	[	A	T	C	G	A	T	C	-	G	A	T	C	G	]
	 * 	qual:	[	>	>	>	>	<	<	<		<	>	>	>	>	]
	 * 	pos:	[	1	2	3		4	5	6	7	8	9	10	11	12	]
	 * 	cigar:	[	M	M	M	I	M	M	M	D	M	M	M	M	M	]
	 * @author covingto
	 *
	 */
	public static class ConformedRead{
		private final SAMRecord rec;
		private final boolean isForward;
		final byte[] ref;
		final byte[] read;
		final byte[] qual;
		final int[] pos;
		final CigarOperator[] cigar;
		private final int mappedpos;
		private final String chr;
		private final int mapEnd;
		private final int mapQual;
		
		public ConformedRead(SAMRecord rec, String chr, int mappedpos, int mapend, byte[] ref, byte[] read, byte[] qual, int[] pos, CigarOperator[] cigar, int mapQual) throws Exception{
			//boolean forward = ! rec.getReadNegativeStrandFlag();
			this(rec, chr, mappedpos, mapend, ref, read, qual, pos, cigar, mapQual, ! rec.getReadNegativeStrandFlag());
		}
		
		public ConformedRead(SAMRecord rec, String chr, int mappedpos, int mapend, byte[] ref, byte[] read, byte[] qual, int[] pos, CigarOperator[] cigar, int mapQual, boolean forward) throws Exception{
			int len = pos.length;
			if (ref.length != len || read.length != len || qual.length != len || cigar.length != len){
				throw new Exception("Arrays are not the same length");
			}
			this.rec = rec;
			this.isForward = forward;
			this.ref = ref;
			this.read = read;
			this.qual = qual;
			this.pos = pos;
			this.cigar = cigar;
			this.chr = chr;
			this.mappedpos = mappedpos;
			this.mapEnd = mapend;
			this.mapQual = mapQual;
		}
		
		
		/**
		 * Generates a new ConformedRead object with all data reflecting the inclusive data of this ConformedRead object where positions are between the indicated start and end ranges.
		 * @param start
		 * @param end
		 * @return
		 * @throws Exception
		 */
		public ConformedRead sliceToPos(int start, int end) throws Exception{
			// find the index for the cut
			if (start < end){ return null; }
			int sstart = 0;
			int send = this.pos.length;
			for (int i = 0; i < this.pos.length; i++){
				if (this.pos[i] == start){
					sstart = i;
					break;
				}
			}
			for (int i = sstart; i < this.pos.length; i++){
				if (this.pos[i] == end){
					send = i;
				} else if (this.pos[i] > end){
					break;
				}
			}
			//return new ConformedRead(this.rec, this.chr, this.pos[sstart], this.pos[send], Arrays.copyOfRange(ref, sstart, send), 
			//		Arrays.copyOfRange(read, sstart, send), Arrays.copyOfRange(qual, sstart, send), Arrays.copyOfRange(pos, sstart, send),
			//		Arrays.copyOfRange(cigar, sstart, send), this.mapQual);
			return sliceToIndex(sstart, send);
					//new ConformedRead(this.chr, this.pos[sstart], this.pos[send], Arrays.copyOfRange(ref, sstart, send), 
					//Arrays.copyOfRange(read, sstart, send), Arrays.copyOfRange(qual, sstart, send), Arrays.copyOfRange(pos, sstart, send),
					//Arrays.copyOfRange(cigar, sstart, send), this.mapQual, this.isForward());
		}
		
		public ConformedRead sliceToIndex(int start, int end) throws Exception{
			return new ConformedRead(this.rec, this.chr, this.pos[start], this.pos[end], Arrays.copyOfRange(ref, start, end), 
					Arrays.copyOfRange(read, start, end), Arrays.copyOfRange(qual, start, end), Arrays.copyOfRange(pos, start, end),
					Arrays.copyOfRange(cigar, start, end), this.mapQual, this.isForward());
		}
		
		public ByteContainer getReadAtGenomicPos(int pos){
			return this.getReadAtGenomicRange(pos, pos);
		}
		
		/**
		 * Return the base quality for this read at (or near) the indicated genomic position.
		 * If there are ins or del after the requested genomic, that quality is returned instead.
		 * @param pos
		 * @return Byte
		 */
		public Byte getQualAtGenomicPos(int pos){
			for (int i = 0; i < this.ref.length; i++){
				if (this.pos[i] == pos){
					if (i < (this.ref.length - 1) && (this.cigar[i+1] == CigarOperator.INSERTION || this.cigar[i+1] == CigarOperator.DELETION)){
						return this.qual[i+1];
					} else {	
						return this.qual[i];
					}
				}
			}
			return null;
		}
		
		public ByteContainer getRefAtGenomicPos(int pos){
			return this.getRefAtGenomicRange(pos, pos);
		}
		
		public int getAlleleEndAtGenomicPos(int end){
			return getAlleleRangeAtGenomicPos(end)[1];
		}
		
		public int getAlleleStartAtGenomicPos(int start){
			return getAlleleRangeAtGenomicPos(start)[0];
		}
		
		/**
		 * ref:			[	a	t	c	g	-	-	-	a	t	c	g	a	t	c	g	]
		 * read:		[	a	t	c	g	a	t	a	a	t	g	c	a	-	-	-	]
		 * pos:			[	1	2	3	4	-1	-1	-1	5	6	7	8	9	10	11	12	]
		 * cigar:		[	m	m	m	m	i	i	i	m	m	m	m	m	d	d	d	]
		 * i:			[	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	]
		 * 
		 * read:		[	a	t	c	t	a	t	a	a	t	g	c	a	-	-	-	]
		 * pos:			[	1	2	3	4	-1	-1	-1	5	6	7	8	9	10	11	12	]
		 * cigar:		[	m	m	m	m	i	i	i	m	m	m	m	m	d	d	d	]
		 * i:			[	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	]
		 * 
		 * read:		[	a	t	c	t	-	-	-	a	t	g	c	a	-	-	-	]
		 * pos:			[	1	2	3	4	-	-	-	5	6	7	8	9	10	11	12	]
		 * cigar:		[	m	m	m	m	-	-	-	m	m	m	m	m	d	d	d	]
		 * i:			[	0	1	2	3	-	-	-	4	5	6	7	8	9	10	11	]
		 * 
		 * say one asks for allele start at position 4.  The answer is subtly different for both reads listed here.
		 * While for the first the allele starts at 4 (it is a match to the reference and is beside an insertion into the reference, 
		 * the results for the second is 3 while the result for the third is again 4.  This is because of the proximity of the nearby indel.
		 * 
		 * This function will return the "most reasonable" / VCF specific genomic position based on the allele present in this read.
		 * @param end
		 * @return
		 */
		public int[] getAlleleRangeAtGenomicPos(int start){
			final int l = this.read.length;
			final int searchSpace = 3;
			int endPosition = -1;
			int startPosition = -1;
			int rightHitOffset = 0;						// records the rightmost position of a non-reference variant for an allele within the search buffer
			int rightHitOffsetl = 0;					// records the leftmost position of an allele to the right of the query position, used for resetting as alleles leave the searchSpace
			int hitOffset = 0;							// records the hit offset of this position lookup
			int refMatchOffset = 0;						// records the leftmost proximal position where the read matches the reference
			// decrement from the end to find the hit offset
			boolean foundHit = false;
			for (int i = l - 1; i > -1; i--){
				// can this be the start of a new allele?
				// for indels the allele always starts with a reference base
				// for snv the allele is the mismatch base
				final int thisPos = this.pos[i];
				final byte ref = this.ref[i];
				final byte read = this.read[i]; 
				final CigarOperator thisCigar = this.cigar[i];
				if (!foundHit && thisPos <= start && thisCigar == CigarOperator.MATCH_OR_MISMATCH){ // case for snp and insertion
					foundHit = true;
					hitOffset = i;
					// ensure that the rightHitOffset and rightHitOffsetl are > hitOffset
					//rightHitOffset = rightHitOffset > hitOffset ? rightHitOffset : hitOffset;
					//rightHitOffsetl = rightHitOffsetl > hitOffset ? rightHitOffsetl : hitOffset;
				} 
				if (! foundHit ){
					// these are the positions to the right of the hit
					if (ref != read){
						rightHitOffsetl = i;
						rightHitOffset = i > rightHitOffset ? i : rightHitOffset; 
					} else {
						// 
						rightHitOffset = i + searchSpace < rightHitOffsetl ? 0 : rightHitOffset;
					}
				} else {
					if (thisCigar == CigarOperator.SOFT_CLIP){
						refMatchOffset = i + 1; // we back up to the last possible hit, this really won't generate an allele
						break; // we can do nothing inside of soft clipping
					}
					// these are the positions to the left of the hit
					if (ref == read) {
						// This looked odd to me too when I was reading through it, but remember that the else clause below resets refMatchOffset to 0
						// in fact we want the maximal position where the ref matches the read so a max here is appropriate
						// also note that we then search against hitOffset, but hitOffset "grows" in the next section (if the if fails)
						refMatchOffset = refMatchOffset > i ? refMatchOffset : i;
						if (i < hitOffset - searchSpace){
							break; // don't need to look any more, save the cycles
						}
					} else {
						// this is interesting, we must be within the buffer to make it this far (we haven't broken yet),
						// but found another indel or variant!! how neat, we get to expand the hit position
						hitOffset = i;
						refMatchOffset = 0; // reset the refmatchoffset 
					}
				}
			}	
			
			// solve the ends of the hits
			if (rightHitOffsetl > hitOffset && rightHitOffset > hitOffset && rightHitOffsetl - hitOffset < searchSpace){
				// indicates that the left end of the right hit is within the search space so we get to set the endPosition
				endPosition = this.pos[rightHitOffset];
			} else {
				// return the position for the hitOffset since this is the first position that we see
				endPosition = this.pos[hitOffset];
			}
			if (rightHitOffset < hitOffset && refMatchOffset + 1 == hitOffset){
				// this is a SNP, we return the actual hit offset as the start position as well
				startPosition = this.pos[hitOffset];
			} else {
				startPosition = this.pos[refMatchOffset];
			}
			//log.log(Level.WARNING, "Conforming positions: refMatchOffset: " + refMatchOffset + " hitOffset: " + hitOffset + " rightHitOffsetl: " + rightHitOffsetl + " rightHitOffset: " + rightHitOffset );
			return new int[] { startPosition, endPosition };
		}
		
		
		/**
		 * Returns a new ByteContainer with the read at the indicated genomic positions.  This will insert .'s (not N's) where the positions would extend beyond the scope of the read.
		 * It may be up to downstream tools to decide how and when to count reads which may have insufficient size to support any particular sequence.
		 * 
		 * Pos:				[	0	1	2	3	4	5	6	7	]
		 * Reference:		[	A	T	C	G	A	T	C	G	]
		 * Read1:			[	A	T	C	A	A	T	C	G	]
		 * Read1:2-6				[	C	A	A	T	C	]
		 * Read2:			[	A	T	C	G	A	T	]
		 * Read2:2-6				[	C	A	A	T	.	]
		 * @param start - genomic start position (inclusive)
		 * @param stop  - genomic end position (inclusive)
		 * @return
		 */
		public ByteContainer getReadAtGenomicRange(int start, int stop){
			// this method returns an array since there might be insertions which do not belong to a genomic position
			boolean hitPos = false;
			List<Byte> bytes = new ArrayList<Byte>((stop - start) + 10);
			for (int i = 0; i < this.read.length; i++){
				// once you hit the position that is always added to the array.  Other switches handle specific cases
				// note that this returns an empty array
				if (this.pos[i] == start){
					hitPos = true;
					bytes.add(this.read[i]);
				} else if (hitPos && this.pos[i] <= stop){ //(hitPos & ( this.cigar[i] == CigarOperator.INSERTION || this.cigar[i] == CigarOperator.DELETION)){
					bytes.add(this.read[i]);
				} else if (this.pos[i] > stop){
					break;
				}
			}
			if (this.pos[this.pos.length - 1] < stop){
				// now fill in the dots
				for (int i = 0; i < (stop - this.pos[this.pos.length - 1]); i++){
					bytes.add(dot);
				}
			}
			return new ByteContainer(bytes.toArray(new Byte[bytes.size()]));
		}
		
		public ByteContainer getRefAtGenomicRange(int start, int stop){
			boolean hitPos = false;
			List<Byte> bytes = new ArrayList<Byte>();
			for (int i = 0; i < this.ref.length; i++){
				// once you hit the position that is always added to the array.  Other switches handle specific cases
				// note that this returns an empty array
				if (this.pos[i] == start){
					hitPos = true;
					bytes.add(this.ref[i]);
				} else if (hitPos && this.pos[i] <= stop){ //(hitPos & ( this.cigar[i] == CigarOperator.INSERTION || this.cigar[i] == CigarOperator.DELETION)){
					bytes.add(this.ref[i]);
				} else if (this.pos[i] > stop){
					break;
				}
			}
			return new ByteContainer(bytes.toArray(new Byte[bytes.size()]));
		}
		
		public CigarOperator[] getCigarAtGenomicRange(int start, int end){
			boolean hitPos = false;
			List<CigarOperator> cigarops = new ArrayList<CigarOperator>((end - start) + 10);
			for (int i = 0; i < this.cigar.length; i++){
				// once you hit the position that is always added to the array.  Other switches handle specific cases
				// note that this returns an empty array
				if (this.pos[i] == start){
					hitPos = true;
					cigarops.add(this.cigar[i]);
				} else if (hitPos && this.pos[i] <= end && this.read[i] != BAMUtils.unk){ //(hitPos & ( this.cigar[i] == CigarOperator.INSERTION || this.cigar[i] == CigarOperator.DELETION)){
					cigarops.add(this.cigar[i]);
				} else if (this.pos[i] > end){
					break;
				}
			}
			return cigarops.toArray(new CigarOperator[cigarops.size()]);
		}
		
		public ByteContainer getQualityAtGenomicRange(int start, int end) {
			boolean hitPos = false;
			List<Byte> bytes = new ArrayList<Byte>((end - start) + 10);
			for (int i = 0; i < this.qual.length; i++){
				// once you hit the position that is always added to the array.  Other switches handle specific cases
				// note that this returns an empty array
				if (this.pos[i] == start){
					hitPos = true;
					bytes.add(this.qual[i]);
				} else if (hitPos && this.pos[i] <= end && this.read[i] != BAMUtils.unk){ //(hitPos & ( this.cigar[i] == CigarOperator.INSERTION || this.cigar[i] == CigarOperator.DELETION)){
					bytes.add(this.qual[i]);
				} else if (this.pos[i] > end){
					break;
				}
			}
			return new ByteContainer(bytes.toArray(new Byte[bytes.size()]));
		}
		
		public Integer[] getINSIndices(){
			List<Integer> ints = new ArrayList<Integer>();
			for (int i = 1; i < this.cigar.length; i++){
				if (this.cigar[i] == CigarOperator.INSERTION & this.cigar[i-1] == CigarOperator.MATCH_OR_MISMATCH){
					ints.add(i-1);
				}
			}
			return ints.toArray(new Integer[ints.size()]);
		}
		
		public Integer[] getDELIndices(){
			List<Integer> ints = new ArrayList<Integer>();
			for (int i = 1; i < this.cigar.length; i++){
				if (this.cigar[i] == CigarOperator.DELETION & this.cigar[i-1] == CigarOperator.MATCH_OR_MISMATCH){
					ints.add(i-1);
				}
			}
			return ints.toArray(new Integer[ints.size()]);
		}
		
		@Override
		public String toString(){
			return this.chr + "\n" + this.mappedpos + "\n" + 
					this.pos[0] + "\t" + this.pos[this.pos.length - 1] + "\n" + new String(this.ref) + "\n" + new String(this.read) + "\n" + 
					StringUtils.join(this.pos, "|");
		}


		public int getMapQuality() {
			return this.mapQual;
		}


		public int readStart() {
			return this.mappedpos;
		}


		public int readEnd() {
			return this.mapEnd;
		}
		
		public ByteContainer getSeqAllele(int start, int end){
			return removePlaceholders(this.getReadAtGenomicRange(start, end));
		}
		
		public static ByteContainer removePlaceholders(ByteContainer seq){
			int redacted = 0;
			for (int i = 0; i < seq.bytes.length; i++){
				if (seq.bytes[i] == unk){
					redacted++;
				}
			}
			if (redacted == 0) { return seq; } // we don't need to redact anything.
			Byte[] result = new Byte[seq.bytes.length - redacted];
			int j = 0;
			for (int i = 0; i < seq.bytes.length; i++){
				if (seq.bytes[i] != unk){
					result[j] = seq.bytes[i];
					j++;
				}
			}
			return new ByteContainer(result);
		}


		public boolean isForward() {
			return this.isForward;
		}

		/**
		 * Check if there is an N within the read for this region
		 * @param start
		 * @param end
		 * @return
		 */
		public boolean containsN(int start, int end){
			ByteContainer refSeq = this.getReadAtGenomicRange(start, end);
			for (int i = 0; i < refSeq.bytes.length; i++){
				if (refSeq.bytes[i] == N){ return true; }
			}
			return false;
		}
		
		/**
		 * Checks if the read and the reference are identical for this genomic range.
		 * @param start
		 * @param end
		 * @return
		 */
		public boolean isReference(int start, int end) {
			ByteContainer refSeq = this.getRefAtGenomicRange(start, end);
			ByteContainer readSeq = this.getReadAtGenomicRange(start, end);
			if (refSeq.bytes.length != readSeq.bytes.length){ return false; }
			for (int i = 0; i < readSeq.bytes.length; i++){
				if (refSeq.bytes[i] != readSeq.bytes[i]){ return false; }
			}
			return true;
		}


		public ByteContainer getSeqQuality(int start, int end) {
			return this.getQualityAtGenomicRange(start, end);
		}


		public String getChr() {
			return this.chr;
		}

		public boolean isRead1() {
			// we could use this.rec.getFirstOfPairFlag() but really don't want to since this shorts if the data aren't really paired
			// instead use 64
			return ((this.rec.getFlags() & 64) != 0); 
		}

	}
	
	public static class Variant{
		public static int MAXMAPQ = 60;
		public class ReadData{
			protected final ByteContainer seqallele;
			protected final ByteContainer refallele;
			protected final byte qual;
			protected final int mapQual;
			protected final int mapPos;
			ReadData(ByteContainer refallele2, ByteContainer seqallele2, byte qual, int mapQual, int mapPos){
				this.refallele = refallele2;
				this.seqallele = seqallele2;
				this.qual = qual;
				this.mapQual = mapQual;
				this.mapPos = mapPos;
			}
			public ByteContainer getSeqAllele(){
				return this.seqallele;
			}
			public ByteContainer getRefAllele(){
				return this.refallele;
			}
			public int getMapPos(){
				return this.mapPos;
			}
		}
		private final String chr;
		private final int pos;
		
		private Set<ByteContainer> variants = new HashSet<ByteContainer>();
		private Set<ByteContainer> references = new HashSet<ByteContainer>();
		private List<ReadData> alleles = new ArrayList<ReadData>();
		
		
		public Variant(String chr, int pos){
			this.chr = chr;
			this.pos = pos;
		}
		
		public void add(ByteContainer refallele, ByteContainer seqallele, byte qual, int mapQual, int mapPos){
			variants.add(seqallele);
			references.add(refallele);
			alleles.add(new ReadData(refallele, seqallele, qual, mapQual, mapPos));
		}

		public String getChr() {
			return chr;
		}

		public int getPos() {
			return pos;
		}

		public Set<ByteContainer> getVariants() {
			return variants;
		}
		
		public Set<ByteContainer> getReferences() {
			return references;
		}

		public int getAlleleCount(String var) {
			int varCount = 0;
			for (ReadData rd : alleles){
				if (rd.seqallele.equals(var)){
					varCount++;
				}
			}
			return varCount;
		}

		public int getTotCount() {
			return alleles.size();
		}
		
		public Byte[] getQuals(){
			Byte[] quals = new Byte[alleles.size()];
			for (int i = 0; i < alleles.size(); i++){
				quals[i] = alleles.get(i).qual;
			}
			return quals;
		}
		
		public Byte[] getAlleleQuals(ByteContainer v){
			List<Byte> quals = new ArrayList<Byte>();
			for (ReadData rd : alleles){
				if (rd.seqallele.equals(v)){
					quals.add(rd.qual);
				}
			}
			return quals.toArray(new Byte[quals.size()]);
		}

		public int getSumAlleleBQ(ByteContainer allele) {
			int sumVarBQ = 0;
			Byte[] quals = this.getAlleleQuals(allele);
			for (int i = 0; i < quals.length; i++){
				sumVarBQ += quals[i].intValue();
			}
			return sumVarBQ;
		}
		
		public Integer[] getMapQuals(){
			Integer[] mapQuals = new Integer[alleles.size()];
			for (int i = 0; i < alleles.size(); i++){
				mapQuals[i] = alleles.get(i).mapQual;
			}
			return mapQuals;
		}
		
		public Integer[] getAlleleMapQuals(Byte[] allele){
			List<Integer> mapQuals = new ArrayList<Integer>();
			for (ReadData rd : alleles){
				if (rd.seqallele.equals(allele)){
					mapQuals.add(rd.mapQual);
				}
			}
			return mapQuals.toArray(new Integer[mapQuals.size()]);
		}
		
		public static String[] getByteListStrs(Collection<Byte[]> sequences){
			List<Byte[]> newsequences = new ArrayList<Byte[]>();
			newsequences.addAll(sequences);
			return getByteListStrs(newsequences);
		}
		
		public static String[] getByteListStrs(List<Byte[]> sequences){
			String[] seqs = new String[sequences.size()];
			for (int i = 0; i < seqs.length; i++){
				Byte[] oldbyte = sequences.get(i);
				byte[] newbyte = new byte[oldbyte.length];
				for (int j = 0; j < oldbyte.length; j++){
					newbyte[j] = oldbyte[j];
				}
				seqs[i] = new String(newbyte);
			}
			return seqs;
		}
		
		@Override
		public String toString(){
			int maxMapQCount = 0;
			int totCount = 0;
			int totQ20Cov = 0;
			for (Integer i : getMapQuals()){
				if (i >= MAXMAPQ){
					maxMapQCount++;
				}
				totCount++;
			}
			for (Byte bq : getQuals()){
				if (bq.intValue() >= 20){
					totQ20Cov++;
				}
			}
			String varCount = "";
			String varQ20Cov = "";
			String sumVarBQ = "";
			for (ByteContainer v : this.getVariants()){
				Byte[] vquals = this.getAlleleQuals(v);
				int thisVarQ20Count = 0;
				int thisVarSumBQ = 0;
				for (Byte bq : vquals){
					thisVarSumBQ += bq.intValue();
					if (bq.intValue() >= 20){
						thisVarQ20Count++;
					}
				}
				varCount += v.toString() + "," + vquals.length + ";";
				varQ20Cov += v.toString() + "," + thisVarQ20Count + ";";
				sumVarBQ += v.toString() + "," + thisVarSumBQ + ";";
			}
			return chr + "\t" + pos + "\t" + StringUtils.join(this.getReferences(), ",") + "\t" + StringUtils.join(this.getVariants(), ",") + "\t" + varCount + 
					"\t" + totCount + "\t" + varQ20Cov + "\t" + totQ20Cov + "\t" + sumVarBQ + 
					"\t" + maxMapQCount;
		}
		
		public static String header(){
			return "chr\tpos\treference\tvariant\tvarCount\totherCount\ttotCount\tvarQ20Cov\ttotQ20Cov\tsumVarBQ\tmaxMapQCount\tvariants";
		}
	}
	
	public static class CRIterator implements Iterator<ConformedRead>, Closeable{
		
		private final String chr;
		private final int start;
		private final int end;
		private final int f;
		private final int F;
		private final IndexedFastaSequenceFile fastaref;
		private final SAMRecordIterator sri;
		
		public CRIterator(SamReader sam, String chr, int start, int end, int f, int F, IndexedFastaSequenceFile fastaref){
			this.sri = sam.query(chr, start, end, false);
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.f = f;
			this.F = F;
			this.fastaref = fastaref;
		}
		
		@Override
		public boolean hasNext() {
			return this.sri.hasNext();
		}
		
		public String getChr() {
		    return this.chr;
		}
		
		public int getStart() {
		    return this.start;
		}
		
		public int getEnd() {
		    return this.end;
		}
		
		@Override
		public ConformedRead next() {
			try{
				final SAMRecord sr = this.sri.next();
				if (((sr.getFlags() & this.F) != 0) // F are the flags to remove: if F and flag are not 0: reject.  F = 0 indicates no filtering
						|| ((sr.getFlags() & this.f) != this.f) // f are the required flags: if f and flag are not f: reject f = 0 indicates no selection
						){
					// log.log(Level.FINEST, "Skipped sam record because of sam flag mask: " + sr.getSAMString());
					return null;
				}
				
				final ConformedRead cr = BAMUtils.conformToReference(sr, fastaref);
				return cr;
			}catch (Exception e){
				log.log(Level.WARNING, "Error proocessing conformed read", e);
				return null;
			}
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void close() throws IOException {
			this.sri.close();
		}
		
	}
	
	public static List<ConformedRead> getConformedReads(BAMInterface bi, String chr, int start, int end, int f, int F, IndexedFastaSequenceFile fastaref){
		final SamReader sam = bi.getSamfilereader();
		List<ConformedRead> reads = getConformedReads(sam, chr, start, end, f, F, fastaref);
		try {
			sam.close();
		} catch (IOException e) {
			log.log(Level.WARNING, "SAM file did not close as expected", e);
		}
		return reads;
	}
	
	public static List<ConformedRead> getConformedReads(SamReader sam, String chr, int start, int end, int f, int F){
		return getConformedReads(sam, chr, start, end, f, F, null);
	}
	
	public static List<ConformedRead> getConformedReads(SamReader sam, String chr, int start, int end, int f, int F, IndexedFastaSequenceFile fastaref){
		List<ConformedRead> reads = new ArrayList<ConformedRead>();
		CRIterator cri = new CRIterator(sam, chr, start, end, f, F, fastaref);
		while (cri.hasNext()){
			try{
				final ConformedRead cr = cri.next();
				if (cr != null){
					reads.add(cr);
				}
			}catch (Exception e){
				log.log(Level.WARNING, "Error proocessing conformed read", e);
				continue;
			}
		}
		try {
			cri.close();
		} catch (IOException e) {
			log.log(Level.SEVERE, "Error in closing conformed read iterator");
			e.printStackTrace();
		}
		return reads;
	}
	
	public static Iterator<ConformedRead> getConformedReadsIterator(SamReader sam, String chr, int start, int end, int f, int F, IndexedFastaSequenceFile fastaref){
		return new CRIterator(sam, chr, start, end, f, F, fastaref);
	}
	
	
	
	public static Variant genotype(SamReader sam, String chr, int pos, String variant) throws Exception{
		return genotype(sam, chr, pos, variant, 0, 1284); /*  1284 = unmapped, not primary align, duplicate*/
	}
	
	public static Variant genotype(SamReader sam, String chr, int pos, String variant, int f, int F) throws Exception{
		Variant var = new Variant(chr, pos);
		List<ConformedRead> reads = BAMUtils.getConformedReads(sam, chr, pos, pos, f, F);
		
		for (ConformedRead cr : reads){
			var.add(cr.getRefAtGenomicPos(pos), cr.getReadAtGenomicPos(pos), cr.getQualAtGenomicPos(pos), cr.getMapQuality(), pos);
		}
		
		return var;
	}
	
	/*
	 * Regexp for MD string.
	 *
	 * \G = end of previous match.
	 * (?:[0-9]+) non-capturing (why non-capturing?) group of digits.  For this number of bases read matches reference.
	 *  - or -
	 * Single reference base for case in which reference differs from read.
	 *  - or -
	 * ^one or more reference bases that are deleted in read.
	 *
	 */
	static final Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");
	

	public static Vector<VariantCall> getVariants(SAMRecord rec){
		return getVariants(rec, null, VarType.ANY);
	}
	
	/*
	 * DEPRICATED
	 */
	public static Vector<VariantCall> getVariants(SAMRecord rec, Integer refpos, VarType type){
		// if refpos is not null then all variants are returned, else we return all reads at the position
		final String md = rec.getStringAttribute(SAMTag.MD.name());
		final Integer chr = rec.getReferenceIndex();
		final int mappedpos = rec.getAlignmentStart();
		final Cigar cigar = rec.getCigar();
		Vector<VariantCall> calls = new Vector<VariantCall>();
		System.out.println("Testing record");
		try{

			int poscount = 0;
			int refcount = 0;

			if (md == null) {
				throw new SAMException("Cannot create reference from SAMRecord with no MD tag, read: " + rec.getReadName());
			}

			byte[] refseq = SequenceUtil.makeReferenceFromAlignment(rec, true);
			byte[] qual = alignQualToReference(rec);
			byte[] newseq = alignSeqToReference(rec);
			if (refseq.length != newseq.length){
				System.out.println("Refseq and Newseq lengths do not match;");
				System.out.println("Refseq len: " + refseq.length + "Newseq len: " + newseq.length);
				System.out.println(refseq);
				System.out.println(newseq);
			}
			for (final CigarElement cigEl : cigar.getCigarElements())
			{
				final int cigElLen = cigEl.getLength();
				final CigarOperator cigElOp = cigEl.getOperator();
				if (cigElOp == CigarOperator.DELETION){
					// build the variant for the deletion
					System.out.println("Found deletion in read " + (mappedpos + refcount));
					int startpos = mappedpos + refcount;
					int endpos = mappedpos + refcount + cigElLen;
					if (refpos == null){
						String reference = new String(Arrays.copyOfRange(refseq, refcount, refcount + cigElLen));
						String variant = "-";
						Byte varQuality = qual[poscount];
						calls.add(new VariantCall(reference, variant, startpos, endpos, chr, varQuality.intValue(), rec.getMappingQuality()));
					} else {
						if (startpos <= refpos & endpos >= refpos & type == VarType.DEL){
							// do a force call and break
							String reference = new String(Arrays.copyOfRange(refseq, refcount, refcount + cigElLen));
							String variant = "-";
							Byte varQuality = qual[poscount];
							calls.add(new VariantCall(reference, variant, startpos, endpos, chr, varQuality.intValue(), rec.getMappingQuality()));
							System.out.println("Added variant to heap");
							break;
						}
					}
					poscount += cigElLen;
					refcount += cigElLen;
				} else if (cigElOp == CigarOperator.SOFT_CLIP){
					// soft clipping will not index up the refcount as these are not counted in the reference map position.
					poscount += cigElLen;
				} else if (cigElOp == CigarOperator.INSERTION){
					// build the variant for the insertion
					// insertions should index up our count of the sequences, but not of the reference map position.
					System.out.println("Found insertion in read " + (mappedpos + refcount));
					int startpos = mappedpos + refcount;
					int endpos = mappedpos + refcount + 1;
					if (refpos == null){
						String variant = new String(Arrays.copyOfRange(newseq, poscount, poscount + cigElLen));
						String reference = "-";
						Byte varQuality = qual[poscount];
						calls.add(new VariantCall(reference, variant, startpos, endpos, chr, varQuality.intValue(), rec.getMappingQuality()));
					} else {
						if (startpos <= refpos & endpos >= refpos & type == VarType.INS){
							String variant = new String(Arrays.copyOfRange(newseq, poscount, poscount + cigElLen));
							String reference = "-";
							Byte varQuality = qual[poscount];
							calls.add(new VariantCall(reference, variant, startpos, endpos, chr, varQuality.intValue(), rec.getMappingQuality()));
							System.out.println("Added variant to heap");
							break;
						}
					}
					poscount += cigElLen;
				} else if (cigElOp == CigarOperator.MATCH_OR_MISMATCH){
					// loop through the array looking for variants
					System.out.println("Found match/mismatch in read " + (mappedpos + refcount));
					if (refpos == null){
					for (int i = 0; i < cigElLen; ++i){
						if (newseq[poscount + i] == refseq[refcount + i]) continue;
						else {
							String reference 	= 	new String(Arrays.copyOfRange(refseq, refcount + i, refcount + i + 1));
							String variant		=	new String(Arrays.copyOfRange(newseq, poscount + i, poscount + i + 1));
							Byte varQuality = qual[poscount + i];
							int startpos 	= 	mappedpos + i;
							int endpos		=	mappedpos + i;
							calls.add(new VariantCall(reference, variant, startpos, endpos, chr, varQuality.intValue(), rec.getMappingQuality()));
						}
					}
					} else {
						int startpos = mappedpos + refcount;
						int endpos = mappedpos + refcount + cigElLen;
						if (startpos <= refpos & endpos >= refpos & type == VarType.SNV){
							int offset = refpos - startpos;
							String reference 	= 	new String(Arrays.copyOfRange(refseq, refcount + offset, refcount + offset + 1));
							String variant		=	new String(Arrays.copyOfRange(newseq, poscount + offset, poscount + offset + 1));
							Byte varQuality = qual[poscount + offset];
							calls.add(new VariantCall(reference, variant, startpos, endpos, chr, varQuality.intValue(), rec.getMappingQuality()));
							break;
						}
					}
					// index up both the ref and the position counts
					poscount += cigElLen;
					refcount += cigElLen;
				}
			} 
		}catch (Exception exc){
			System.out.println("Caught exception while processing read.");
			System.out.println(exc.toString());
		}
		return calls;
	}
	
	public static ConformedRead conformToReference(SAMRecord rec) throws Exception{
		return conformToReference(rec, null);
	}
	
	/**
	 * 
	 * Results in a ConformedRead object, who's main elements are sets of byte[] and enum[] arrays that line up the read with the reference.
	 * <pre>
	 * Ex;
	 * 	read;	ATCGATCGATCG 
	 * 	qual;	>>>><<<<>>>>	
	 * 	CIGAR;	3M1I3M1D5M
	 * 	pos;	1
	 * to 
	 * 	ref:	[	A	T	C	-	A	T	C	?	G	A	T	C	G	]
	 * 	read:	[	A	T	C	G	A	T	C		G	A	A	C	G	]
	 * 	qual:	[	>	>	>	>	<	<	<		<	>	>	>	>	]
	 * 	pos:	[	1	2	3	3	4	5	6	7	8	9	10	11	12	]
	 * 	cigar:	[	M	M	M	I	M	M	M	D	M	M	M	M	M	]
	 * </pre>
	 * @param rec
	 * @return
	 * @throws Exception 
	 */
	public static ConformedRead conformToReference(SAMRecord rec, IndexedFastaSequenceFile fastaref) throws Exception{
		final byte[] seq = rec.getReadBases();
		final byte[] qual = rec.getBaseQualities();
		// we must fill in deletions on our own
		final Cigar cigar = rec.getCigar();
		int maxOutputLength = 0;
		final int mapPos = rec.getUnclippedStart();
		if (cigar == null) {
			throw new SAMException("Cannot create reference from SAMRecord with no CIGAR, read: " + rec.getReadName());
		}
		for (final CigarElement cigarElement : cigar.getCigarElements()) {
			final CigarOperator cigElOp = cigarElement.getOperator();
			if (cigElOp == CigarOperator.HARD_CLIP || cigElOp == CigarOperator.PADDING){ continue; }
			maxOutputLength += cigarElement.getLength();
		}
		byte[] refseq;
		boolean fromRead;
		if (fastaref != null){
			fromRead = false;
			refseq = SynchronousIndexedFastaReader.getSubsequenceAt(fastaref, rec.getReferenceName(), mapPos, mapPos + maxOutputLength).getBases();
		} else {
			refseq = SequenceUtil.makeReferenceFromAlignment(rec, true);
			fromRead = true;
			// note that the refseq array contains '-' where there are insertions and '0' where there is soft clipping
		}
		
		// make the conformed sequences
		byte[] refArray = new byte[maxOutputLength];
		byte[] readArray = new byte[maxOutputLength];
		byte[] qualArray = new byte[maxOutputLength];
		int[] posArray = new int[maxOutputLength];
		CigarOperator[] cigarArray = new CigarOperator[maxOutputLength];
		
		// trackers
		int seqPos = 0; // this tracker will maintain our index in the sequence and the qualities.  This indexes at insertions but not at deletions.
		int refPos = 0; // this tracker will maintain the index of the reference.  This indexes as deletions but not at insertions.
		int arrayPos = 0; // this tracker indicates where we are writing in the arrays.  This indexes at everything.
		int gpos = mapPos;
		for (final CigarElement cigEl : cigar.getCigarElements()){
			final int cigElLen = cigEl.getLength();
			final CigarOperator cigElOp = cigEl.getOperator();
			try {
				if (cigElOp == CigarOperator.DELETION){
					// deletion from the reference, pos data is consistent with the maped pos
					for (int i = 0; i < cigElLen; i++){
						refArray[arrayPos + i] = refseq[refPos + i];
						qualArray[arrayPos + i] = unk;
						readArray[arrayPos + i] = unk;
						posArray[arrayPos + i] = gpos;
						gpos++;
						cigarArray[arrayPos + i] = cigElOp;
					}
					arrayPos += cigElLen;
					refPos += cigElLen;
				} else if (cigElOp == CigarOperator.N){
					// skips (introns according to the spec) aren't really deletions in the genomic sense so calling them unk in our vocabulary is not correct.
					// instead they are labeled dot so that other applications can deal with them.
					for (int i = 0; i < cigElLen; i++){
						refArray[arrayPos + i] = refseq[refPos + i];
						qualArray[arrayPos + i] = dot;
						readArray[arrayPos + i] = dot;
						posArray[arrayPos + i] = gpos;
						gpos++;
						cigarArray[arrayPos + i] = cigElOp;
					}
					arrayPos += cigElLen;
					refPos += cigElLen;
				} else if (cigElOp == CigarOperator.INSERTION){
					for (int i = 0; i < cigElLen; i++){
						posArray[arrayPos + i] = gpos;
						if (fromRead) { refArray[arrayPos + i] = refseq[refPos + i]; }
						else { refArray[arrayPos + i] = unk; }
						qualArray[arrayPos + i] = qual[seqPos + i];
						readArray[arrayPos + i] = seq[seqPos + i];
						cigarArray[arrayPos + i] = cigElOp;
					}
					if (fromRead){ refPos += cigElLen; }
					arrayPos += cigElLen;
					seqPos += cigElLen;
				} else if (cigElOp == CigarOperator.MATCH_OR_MISMATCH || cigElOp == CigarOperator.SOFT_CLIP || cigElOp == CigarOperator.EQ || cigElOp == CigarOperator.X){
					for (int i = 0; i < cigElLen; i++){
						// everything updates
						refArray[arrayPos + i] = refseq[refPos + i];
						qualArray[arrayPos + i] = qual[seqPos + i];
						readArray[arrayPos + i] = seq[seqPos + i];
						posArray[arrayPos + i] = gpos;
						gpos++;
						cigarArray[arrayPos + i] = cigElOp;
					}
					refPos += cigElLen;
					arrayPos += cigElLen;
					seqPos += cigElLen;
				} else if (cigElOp == CigarOperator.HARD_CLIP || cigElOp == CigarOperator.PADDING) { 
					continue; 
				} else {
					log.log(Level.SEVERE, "Found unknown operator " + cigElOp);
					throw new Exception("Unknown operator in read: " + rec.toString());
				}
			} catch (Exception e) {
				log.log(Level.WARNING, "Caught exception parsing conformed read.");
				log.log(Level.WARNING, "SEQ:  " + new String(seq));
				log.log(Level.WARNING, "RSEQ: " + new String(refseq));
				log.log(Level.WARNING, "RefArray:  " + new String(refArray));
				log.log(Level.WARNING, "ReadArray: " + new String(readArray));
				log.log(Level.WARNING, "CigarOp: " + cigElOp.toString() + " Length: " + cigElLen);
				throw e;
			}
		}
		return new ConformedRead(rec, rec.getReferenceName(), mapPos, rec.getUnclippedEnd(), refArray, readArray, qualArray, posArray, cigarArray, rec.getMappingQuality());
	}
	
	public static byte[] alignToReference(byte[] seq, Cigar cigar, SAMRecord rec){
		int maxOutputLength = 0;
		int poscount = 0;
		int seqcount = 0;

		if (cigar == null) {
			throw new SAMException("Cannot create reference from SAMRecord with no CIGAR, read: " + rec.getReadName());
		}
		for (final CigarElement cigarElement : cigar.getCigarElements()) {
			maxOutputLength += cigarElement.getLength();
		}
		byte[] newseq = new byte[maxOutputLength];

		// move through the cigar string to line up the bases, insertions are considered variants out of the box so those can be processed here.
		for (final CigarElement cigEl : cigar.getCigarElements())
		{
			final int cigElLen = cigEl.getLength();
			final CigarOperator cigElOp = cigEl.getOperator();

			if (cigElOp == CigarOperator.SKIPPED_REGION) {
				throw new SAMException("Found skipped regions in " + rec.getReadName() + " cigar string; " + cigar);
			} else if (cigElOp == CigarOperator.DELETION){
				for (int i = 0; i < cigElLen; ++i){
					newseq[poscount++] = N;
				}
			} else {
				for (int i = 0; i < cigElLen; ++i){
					newseq[poscount++] = seq[seqcount++];
				}
			}

		}
		return newseq;
	}

	public static byte[] alignSeqToReference(SAMRecord rec){
		final Cigar cigar = rec.getCigar();
		final byte[] seq = rec.getReadBases();
		return alignToReference(seq, cigar, rec);
	}
	
	public static byte[] alignQualToReference(SAMRecord rec){
		final Cigar cigar = rec.getCigar();
		final byte[] qual = rec.getBaseQualities();
		return alignToReference(qual, cigar, rec);
	}
	
	public static String byteArrayToString(Byte[] array){
		byte[] newBytes = new byte[array.length];
		for (int i = 0; i < array.length; i++){
			newBytes[i] = array[i].byteValue();
		}
		return new String(newBytes);
	}
	
	public static void main(String[] args) throws Exception{
		// testing
		System.out.println("Testing BAMUtils.  This is for debugging only!");
		String inBamName = args[0];
		String chr = args[1];
		Integer start = Integer.parseInt(args[2]);
		Integer end = Integer.parseInt(args[3]);
		String variant = args[4];
		//VarType vartype = VarType.valueOf(args[5]);
		
		SamReader samReader = SamReaderFactory.makeDefault().open(new File(inBamName));
		List<ConformedRead> conformedReads = getConformedReads(samReader, chr, start, end, 0, 0);
		for (ConformedRead c : conformedReads){
			System.out.println("NewRead");
			System.out.println("Reference:    " + new String(c.ref));
			System.out.println("Read:         " + new String(c.read));
			byte[] qualcopy = Arrays.copyOf(c.qual, c.qual.length);
			for (int i = 0; i < qualcopy.length; i++){
				qualcopy[i] = (byte) (qualcopy[i] + 33);
			}
			System.out.println("Qual:         " + new String(qualcopy));
		}
		Variant testAllele = genotype(samReader, chr, start, variant);
		System.out.println(Variant.header());
		System.out.println(testAllele.toString());
		
		
	}
	
}

