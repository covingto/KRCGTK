package org.bcm.hgsc.cancer.pacbio;

import org.bcm.hgsc.cancer.utils.Orientation;

public class ReferenceSequence {
	private final String chr;
	private final Integer start;
	private final Integer end;
	private final String s;
	private final Orientation o;
	
	public ReferenceSequence(String chr, Integer start, Integer end, String seq, Orientation o){
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.s = seq;
		this.o = o;
	}
	
	public String name(){
		return this.chr + "_" + this.start + "_" + this.end + "_" + this.o;
	}
	
	public ReferenceSequence subsequence(Integer start, Integer end) throws Exception{
		String newRef;
		Integer nstart;
		Integer nend;
		
		switch (this.o){
			case POS:
				newRef = (String) this.s.subSequence(start, end);
				nstart = this.start + start;
				nend = this.end - (this.s.length() - end);
				break;
			case NEG:
				newRef = (String) this.s.subSequence(start, end);
				nstart = this.start + (this.s.length() - end);
				nend = this.end - start;
				break;
			default:
				throw new Exception("Incompatible argument, Orientation neither positive nor negative");
		}
		
		return new ReferenceSequence(this.chr, nstart, nend, newRef, this.o);
	}
}
