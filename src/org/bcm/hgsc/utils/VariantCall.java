package org.bcm.hgsc.utils;

public class VariantCall {
	private final String reference;
	private final String variant;
	private final int startpos;
	private final int endpos;
	private final int chrom;
	private final int baseQuality;
	private final int mapQuality;
	
	public VariantCall(String ref, String var, int start, int end, int chr, int baseQuality, int mapQuality){
		this.reference = ref;
		this.variant = var;
		this.startpos = start;
		this.endpos = end;
		this.chrom = chr;
		this.baseQuality = baseQuality;
		this.mapQuality = mapQuality;
	}
	
	public int getBaseQuality() {
		return this.baseQuality;
	}
	
	public int getMapQuality() {
		return this.mapQuality;
	}

	public String getReference() {
		return reference;
	}

	public String getVariant() {
		return variant;
	}

	public int getStartpos() {
		return startpos;
	}

	public int getEndpos() {
		return endpos;
	}

	public int getChrom() {
		return chrom;
	}
	
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + chrom;
		result = prime * result + endpos;
		result = prime * result
				+ ((reference == null) ? 0 : reference.hashCode());
		result = prime * result + startpos;
		result = prime * result + ((variant == null) ? 0 : variant.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		//System.out.println("Running Equals");
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		VariantCall other = (VariantCall) obj;
		if (chrom != other.chrom)
			return false;
		if (endpos != other.endpos)
			return false;
		if (reference == null) {
			if (other.reference != null)
				return false;
		} else if (!reference.equals(other.reference))
			return false;
		if (startpos != other.startpos)
			return false;
		if (variant == null) {
			if (other.variant != null)
				return false;
		} else if (!variant.equals(other.variant))
			return false;
		return true;
	}
	
	public String toString(){
		return "Chrom; " + this.chrom + " Start; " + this.startpos + " End; " + this.endpos + " Ref; " + this.reference + " Variant; " + this.variant;
	}
	
}
