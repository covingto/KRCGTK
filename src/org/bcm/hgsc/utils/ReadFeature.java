package org.bcm.hgsc.utils;

import org.bcm.hgsc.cancer.bed.BEDRegion;
import org.bcm.hgsc.cancer.utils.Chromosome;

public class ReadFeature extends BEDRegion{
	public final Chromosome 	chr;
	public final String			read;
	
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + ((chr == null) ? 0 : chr.hashCode());
		result = prime * result + ((read == null) ? 0 : read.hashCode());
		result = prime * result + start;
		result = prime * result + stop;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (!(obj instanceof ReadFeature))
			return false;
		ReadFeature other = (ReadFeature) obj;
		if (chr == null) {
			if (other.chr != null)
				return false;
		} else if (!chr.equals(other.chr))
			return false;
		if (read == null) {
			if (other.read != null)
				return false;
		} else if (!read.equals(other.read))
			return false;
		if (start != other.start)
			return false;
		if (stop != other.stop)
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "ReadFeature [read=" + read + ", chrom=" + chrom + ", start="
				+ start + ", stop=" + stop + "]";
	}

	public ReadFeature(Chromosome chr, Integer p, String r){
		super(chr.toString(), p, p);
		this.chr = chr;
		this.read = r;
	}

	public int compareTo(ReadFeature br) {
		if (this.chrom.equals(br.chrom)){
			if (this.start == br.start){
				return this.read.compareTo(br.read);
			} else {
				return this.start - br.start;
			}
		} else {
			return this.chrom.compareTo(br.chrom);
		}
	}
	
	public boolean near(ReadFeature rf, int buffer){
		if (this.chrom.equals(rf.chrom)){
			if (Math.abs(this.start - rf.start) < buffer){
				return true;
			}
		}
		return false;
	}
	
	
}
