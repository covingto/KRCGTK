package org.bcm.hgsc.cancer.utils;


public class Chromosome {
	private final Integer chrnum;
	private final String 	chrorg;
	public Chromosome(String chr){
		String chrnew = chr.replace("chr", "");
		Integer i = 0;
		try {
			i = Integer.parseInt(chrnew);
		} catch (NumberFormatException e) {
			if (chrnew.equals("X")){
				i = 23;
			} else 
			if (chrnew.equals("Y")) {
				i = 24;
			} // else
			//if (Settings.debug){
			//	System.out.println("Error in decoding chr: " + chr + ", set to " + i);
			//	}
		}
		this.chrnum = i;
		this.chrorg = chr;
	}
	
	public boolean isGreaterThan(Chromosome test){
		return this.chrnum > test.chrnum;
	}
	
	public boolean equals(Chromosome test){
		return this.chrnum == test.chrnum;
	}
	
	public int compare(Chromosome test){
		return this.chrnum - test.chrnum;
	}
	
	public String toString(){
		return this.chrorg;
	}
	
	public int getChrnum(){
		return this.chrnum;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chrnum == null) ? 0 : chrnum.hashCode());
		result = prime * result + ((chrorg == null) ? 0 : chrorg.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof Chromosome))
			return false;
		Chromosome other = (Chromosome) obj;
		if (chrnum == null) {
			if (other.chrnum != null)
				return false;
		} else if (!chrnum.equals(other.chrnum))
			return false;
		if (chrorg == null) {
			if (other.chrorg != null)
				return false;
		} else if (!chrorg.equals(other.chrorg))
			return false;
		return true;
	}
	
	
}
