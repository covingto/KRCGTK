package org.bcm.hgsc.cancer.bed;

public class BEDRegion implements Comparable<BEDRegion> {
	protected final String chrom;
	protected final int start;
	protected final int stop;
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chrom == null) ? 0 : chrom.hashCode());
		result = prime * result + start;
		result = prime * result + stop;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof BEDRegion))
			return false;
		BEDRegion other = (BEDRegion) obj;
		if (chrom == null) {
			if (other.chrom != null)
				return false;
		} else if (!chrom.equals(other.chrom))
			return false;
		if (start != other.start)
			return false;
		if (stop != other.stop)
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "BEDRegion [chrom=" + chrom + ", start=" + start + ", stop="
				+ stop + "]";
	}

	public BEDRegion(String chrom, int i, int j) {
		// TODO Auto-generated constructor stub
		this.chrom = chrom;
		this.start = i;
		this.stop = j;
	}

	public String getSequence() {
		// TODO Auto-generated method stub
		return this.chrom;
	}

	public int getStart() {
		// TODO Auto-generated method stub
		return this.start;
	}

	public int getEnd() {
		// TODO Auto-generated method stub
		return this.stop;
	}

	@Override
	public int compareTo(BEDRegion br) {
		if (this.chrom.equals(br.chrom)){
			if (this.start == br.start){
				return this.stop - br.stop;
			} else {
				return this.start - br.start;
			}
		} else {
			return this.chrom.compareTo(br.chrom);
		}
	}

}
