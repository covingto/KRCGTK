package org.bcm.hgsc.cancer.pacbio;

import java.util.ArrayList;
import java.util.List;

import org.bcm.hgsc.cancer.utils.Orientation;

/**
 * Synthetic references are collections of reference sequences organized into one sequence.
 * @author covingto
 *
 */
public class SyntheticReference {
	//private final String name;
	private final String s;
	private final Integer clipbuffer; // limit beyond which clipping is not allowed
	private final List<ReferencePosition> positions;
	private boolean cantclip = false;
	private final Integer clipstart; // the index in s for where to start the clipping, supplied at creation time 
	//private final String chr;
	//private final Integer start;
	//private final Integer end;
	//private final Orientation o;

	private class ReferencePosition {
		private final String chr;
		private final Integer start;
		private final Integer end;
		private final Orientation o;
		
		public ReferencePosition(String chr, Integer start, Integer end, Orientation o){;
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.o = o;
		}
	}
	
	/**
	 * Generates a blank or empty SyntheticReference, this reference can not be clipped.
	 */
	public SyntheticReference(){
		this(new ArrayList<ReferencePosition>(), "");
	}
	
	public SyntheticReference(List<ReferencePosition> positions, String s){
		this(positions, s, 0, 0);
	}
	
	public SyntheticReference(List<ReferencePosition> positions, String s, Integer clipbuffer, Integer clipstart){
		this.s = s;
		this.clipbuffer = clipbuffer;
		this.positions = positions;
		this.clipstart = clipstart;
	}
	
	public void setCantClip(){
		this.cantclip = true;
	}
	
	public SyntheticReference clip(Integer leftclip) throws ExcessiveClipException, DoubleClipException {
		// return a new synthetic referece that is clipped in on the right hand side.
		if (-leftclip > this.clipbuffer){
			throw new ExcessiveClipException("Clip length exception: " + -leftclip + " > " + this.clipbuffer);
		}
		if (this.cantclip) {
			throw new DoubleClipException("This reference as already been clipped");
		}
		// copy the positions, except the last
		List<ReferencePosition> pos = new ArrayList<ReferencePosition>();
		for (int i = 0; i < positions.size() - 1; i++){
			pos.add(positions.get(i));
		}
		ReferencePosition lastpos = positions.get(positions.size() - 1);
		
		/** set coordinates for the new reference position, this should adjust the positions to reflect the start and end of the genomic sequence that is included 
		 *  in the new synthetic reference.
		 *  
		 *  ===========================================================================================
		 *  |                      |                    |                                             |
		 *  lastpos.start		   sstart				send										  lastpos.end
		 *  
		 *  											|---------------------------------------------|
		 *  														s.length() - this.clipstart()
		 */   
		ReferencePosition newlastpos;
		if (lastpos.o == Orientation.POS){
			newlastpos = new ReferencePosition(lastpos.chr, lastpos.start, lastpos.end - (this.s.length() - this.clipstart) + leftclip, lastpos.o);
		} else {
			newlastpos = new ReferencePosition(lastpos.chr, lastpos.start + (this.s.length() - this.clipstart) - leftclip, lastpos.end, lastpos.o);
		}
		pos.add(newlastpos);
		
		System.out.println("Clipping from " + this.clipstart + " by " + leftclip);
		
		System.out.println("Clipped sequences from :");
		System.out.println(poslistString(this.positions));
		System.out.println("To:");
		System.out.println(poslistString(pos));
		
		
		// add the new position to the growing list of positions
		
		// generate the new string, this is actually easier than it looks
		String newrefstring = (String) this.s.subSequence(0, this.clipstart + leftclip);
		SyntheticReference newreference = new SyntheticReference(pos, newrefstring);
		newreference.setCantClip();
		return newreference;
	}
	
	private static String poslistString(List<ReferencePosition> pos){
		String name = "";
		for (ReferencePosition rp : pos) {
			name = name + ">" + rp.chr + "_" + rp.start + "_" + rp.end + "_" + rp.o; 
		}
		return name;
	}
	
	/**
	 * Returns a synthetic reference object representing the addition of the new reference position and sequence to the SyntheticReference object.
	 * This function also sets the clipbuffer (maximum amount of sequence that can be clipped) and the clipposition (location where clipping begins).
	 * @param rp
	 * @param seq
	 * @param clipbuffer
	 * @param seqbuffer
	 * @return SyntheticReference - newsynthref
	 * @throws DoubleClipException
	 */
	public SyntheticReference add(ReferencePosition rp, String seq, Integer clipbuffer, Integer seqbuffer) throws DoubleClipException{
		// return a new synthetic reference that has a genomic position added to it
		
		// build the position array
		if (!this.cantclip){ throw new DoubleClipException("Adding to synthetic reference that can still be clipped is not allowed."); }
		List<ReferencePosition> pos = this.positions;
		pos.add(rp); // sure throwing these things around might seem like a bad idea, but all of the rp attributes are final so that will save us from modification.
		
		// build the new string
		String newrefstring = this.s + seq;
		System.out.println("Added new references to the sequence");
		// System.out.println(this.s + "<||||||>" + seq);
		System.out.println("Clipbuffer: " + clipbuffer + ", Seqbuffer: " + seqbuffer);
		System.out.println("Setting clipposition to: " + (newrefstring.length() - seqbuffer));
		SyntheticReference newref = new SyntheticReference(pos, newrefstring, clipbuffer, newrefstring.length() - seqbuffer);
		return newref;
	}

	/**
	 * Returns the name that references this synthetic reference, this is already in fasta format.  The name is the concatenation of the chr, start, end, and orientations of all reference positions used for this synthetic reference.
	 * @return String - name
	 */
	public String fastaName() {
		String name = "";
		for (ReferencePosition rp : positions) {
			name = name + ">" + rp.chr + "_" + rp.start + "_" + rp.end + "_" + rp.o; 
		}
		return name;
	}

	public String fastaSeq() {
		return this.s;
	}

	/**
	 * Adds a new sequence to the growing synthetic reference.
	 * Left clip is the amount of sequence that we need to move to the left (5') end of the reference sequence, a positive number results in reducing the genomic position.
	 * Right clip is the amount of sequence that we need to move to the right (3') end of the reference sequence, a positive number results in increasing the genomic position.
	 * This sets the clipbuffer (internal parameter controlling amount of clipping of the sequence) to 50% of the reference sequence length matched by the query.
	 * This sets the clipstart to the length of the combined query minus the right clip size.
	 * @param br
	 * @param leftclip
	 * @param rightclip
	 * @throws Exception
	 */
	public SyntheticReference add(BlastRow br, Integer leftclip, Integer rightclip) throws Exception {
		Integer rpstart = br.o == Orientation.POS ? br.sstart - leftclip : br.send + leftclip;
		Integer rpend	= br.o == Orientation.POS ? br.send + rightclip : br.sstart - rightclip;
		ReferencePosition rp = new ReferencePosition(br.schr, rpstart, rpend, br.o);
		
		// get the new sequence
		String newseq = br.getSequence(leftclip, rightclip);
		// set the clipbuffer to 50% of insertion section of the read
		Integer clipbuffer = Math.abs((br.send - br.sstart)/2);
		return this.add(rp, newseq, clipbuffer, rightclip);
	}
}
