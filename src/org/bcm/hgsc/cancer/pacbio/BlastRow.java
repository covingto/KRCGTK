package org.bcm.hgsc.cancer.pacbio;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;

import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.cancer.utils.Chromosome;
import org.bcm.hgsc.cancer.utils.Orientation;

public class BlastRow {
	private final String string;
	public enum Unique { AMBIGUOUS, UNIQUE, REDUNDANT };
	// q (query) indicates the read
	public final Integer qstart;
	public final Integer qend;
	// s (subject) indicates the reference
	public final String schr;
	public final Integer sstart;
	public final Integer send;
	public final Orientation o;
	public final String q;
	public static final Unique AMBIGUOUS = Unique.AMBIGUOUS;
	public static final Unique UNIQUE = Unique.UNIQUE;
	public static final Unique REDUNDANT = Unique.REDUNDANT;
	public final IndexedFastaSequenceFile reffasta;
	
	public Unique unique = Unique.UNIQUE;
	
	public BlastRow(String s, IndexedFastaSequenceFile reffasta){
		this(s.split("\t"), reffasta);
	}
	
	public BlastRow(String q, Integer qstart, Integer qend, String schr, Integer sstart, Integer send, IndexedFastaSequenceFile reffasta, String str){
		this.q = q;
		this.reffasta = reffasta;
		this.qstart = qstart;
		this.qend = qend;
		this.schr = schr;
		this.sstart = sstart;
		this.send = send;
		if (this.send > this.sstart){
			this.o = Orientation.POS;
		} else {
			this.o = Orientation.NEG;
		}
		this.string = str;
	}
	
	public BlastRow(String[] vals, IndexedFastaSequenceFile reffasta) throws NumberFormatException, IndexOutOfBoundsException{
		this(vals[0], 
				Integer.parseInt(vals[6]), // query start position
				Integer.parseInt(vals[7]), // query end position
				vals[1], // chromosome match
				Integer.parseInt(vals[8]), // subject (reference) start
				Integer.parseInt(vals[9]), // subject (reference) end
				reffasta, 
				StringUtils.join(vals, "\t"));
	}
	
	@Override
	public String toString(){
		return this.string;
	}
	
	public BlastRow invert(){
		return new BlastRow(this.q, this.qend, this.qstart, this.schr, this.send, this.sstart, this.reffasta, this.string);
	}
	
	public String getSequence(Integer leftclip, Integer rightclip) throws Exception{
		// get the sequence relative to the read position
		
		switch (this.o){
			case POS:
				try {
					return getGenomicSequence(this.schr, this.sstart - leftclip, this.send + rightclip).toString();
				} catch (Exception e) {
					System.out.println("Sequence exception " + "Map start: " + this.sstart + ", Map End: " + this.send + 
							", Leftclip: " + leftclip + ", Rightclip: " + rightclip + ", Orientation: " + this.o);
					throw e;
				}
			case NEG:
				try {
					return SequenceUtil.reverseComplement(getGenomicSequence(this.schr, this.send - rightclip, this.sstart + leftclip).toString());
				} catch (Exception e) {
					System.out.println("Sequence exception " + "Map start: " + this.sstart + ", Map End: " + this.send + 
							", Leftclip: " + leftclip + ", Rightclip: " + rightclip + ", Orientation: " + this.o);
					throw e;
				}
			default:
				throw new Exception ("Orientation neither POS nor NEG");
		}
	}
	
	private String getGenomicSequence(String chr, Integer start, Integer end){
		if ( this.reffasta == null ){ return null; }
		return new String(this.reffasta.getSubsequenceAt(chr, start, end).getBases());
	}

	public boolean isGreaterThan(BlastRow rb) {
		// compare the chromosomes first
		Chromosome ca = new Chromosome(this.schr);
		Chromosome cb = new Chromosome(rb.schr);
		int diff = ca.compare(cb);
		if (diff > 0) {
			// ca is greater in coordinate than cb
			return true;
		} else if (diff < 0){
			// ca is less than cb
			return false;
		} else {
			int posdiff = ((this.sstart + this.send) / 2 ) - ((rb.sstart + rb.send) / 2);
			if (posdiff > 0) {
				return true;
			} else {
				return false;
			}
		}
	}
	
	public static Unique unique(BlastRow m1, BlastRow m2){
		Float m1ol = BlastRow.overlapFrc(m1.qstart, m1.qend, m2.qstart, m2.qend);
		Float m2ol = BlastRow.overlapFrc(m2.qstart, m2.qend, m1.qstart, m1.qend);
		if ( m1ol > 0.90 & m2ol > 0.90 ) { 
			return BlastRow.AMBIGUOUS;
		} 
		return BlastRow.UNIQUE;
	}
	
	/**
	 * Calculate the overlap of the reads.
	 * READ		------------------------------------------------------------------------------
	 * MAP1				s1---------------------------------------e1
	 * MAP2				s2---------------------------------------e2 100% overlap with MAP1, both should be ambiguous
	 * MAP3													s2-----------------------------e2 < 90% overlap, considered unique
	 * MAP4						s2-------e2							entirely contained within another read, considered redundant
	 * 
	 * @param s1
	 * @param e1
	 * @param s2
	 * @param e2
	 * @return
	 */
	public static Float overlapFrc(Integer s1, Integer e1, Integer s2, Integer e2){
		// calculate the size of the first mapping
		Integer len = e1 - s1;
		
		// overlap area is the size of the bases that are contained
		Integer ol = Math.min(e1,  e2) - Math.max(s1,  s2);
		return Float.valueOf(ol) / Float.valueOf(len);
	}
}
