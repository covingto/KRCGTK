package org.bcm.hgsc.cancer.sv;

import java.util.Comparator;

import org.bcm.hgsc.cancer.pacbio.BlastRow;
import org.bcm.hgsc.cancer.utils.Chromosome;
import org.bcm.hgsc.cancer.utils.Orientation;
import org.bcm.hgsc.utils.ReadFeature;

/**
 * Junction class is a container for junctions in genomes.
 * @author covingto
 *
 */
public class Junction implements Comparable<Junction> {
	public static enum FusionType{ INTERCHROM, INTRACHROM };
	public final ReadFeature a;
	public final ReadFeature b;
	public final Orientation ol;
	public final Orientation or;
	// the blast entry for the genome normalized right side of the junction
	public final BlastRow brr;
	// the blast entry for the genome normalized left side of the junction 
	public final BlastRow brl;
	public final static Integer buffer = 500;
	public final String read;
	public final String zmw;
	public final FusionType fusionType;
	// if inverted read space info should reverse brl and brr
	public final boolean invert;
	// compareables
	public static Comparator<Junction> jcomparable = new Comparator<Junction>() {

		@Override
		public int compare(Junction arg0, Junction arg1) {
			return arg0.compareTo(arg1);
		}

	};
	
	public static Comparator<Junction> qstartcompare = new Comparator<Junction>() {
		
		@Override
		public int compare(Junction arg0, Junction arg1) {
			int minj1 = Math.min(Math.min(arg0.brl.qstart, arg0.brl.qend), Math.min(arg0.brr.qstart, arg0.brr.qend));
			int minj2 = Math.min(Math.min(arg1.brl.qstart, arg1.brl.qend), Math.min(arg1.brr.qstart, arg1.brr.qend));
			return minj1 - minj2;
			}
		};
	
	public Junction(Chromosome chrl, Chromosome chrr, Integer posl, Integer posr, Orientation ol, Orientation or, BlastRow brl, BlastRow brr, boolean invert){
		this.a = new ReadFeature(chrl, posl, brl.q);
		this.b = new ReadFeature(chrr, posr, brr.q);
		this.ol = ol;
		this.or = or;
		this.read = brl.q;
		this.zmw = readToZMW(this.read);
		this.brl = brl;
		this.brr = brr;
		if (chrl.equals(chrr)){
			this.fusionType = FusionType.INTRACHROM;
		} else {
			this.fusionType = FusionType.INTERCHROM;
		}
		this.invert = invert;
	}
	
	public static String readToZMW(String s){
		String[] zmwl = s.split("/");
		return zmwl[0].substring(zmwl[0].length() - 20, zmwl[0].length()) + "/" + zmwl[1];
	}
	
	public double gapRatio(){
		if (this.fusionType.equals(FusionType.INTERCHROM)){
			return Float.MAX_VALUE;
		} else {
			double readspace = this.brr.qstart - this.brl.qend;
			double seqspace  = this.brr.sstart - this.brl.send;
			return seqspace / readspace;
		}
	}
	
	/**
	 * returns a long for the gap (extra or removed in some cases) bases in the read.
	 * @return
	 */
	public long readGap(){
		return this.queryPosRight() - this.queryPosLeft();
	}
	
	/** 
	 * returns a long for the space in the reference.  
	 * @return
	 */
	public long referenceGap(){
		if (this.fusionType.equals(FusionType.INTERCHROM)){
			return Integer.MAX_VALUE;
		} else {
			return this.brr.sstart - this.brl.send;
		}
	}
	
	public boolean like(Junction test, Integer buffer){
		return this.a.near(test.a, buffer) && this.b.near(test.b, buffer);
		
	}

	public String juctionOuptut() {
		return this.a.chr + "\t" + this.a.getStart() + "\t" + this.ol + "\t" + this.b.chr + "\t" + this.b.getStart() + "\t" + this.or + "\t" + this.a.read + "\t" + this.zmw + "\t" + this.invert;
	}

	@Override
	public int compareTo(Junction o) {
		int compa = this.a.compareTo(o.a);
		int compb = this.b.compareTo(o.b);
		if (compa == 0){
			return compb;
		} else {
			return compa;
		}
	}
	
	/**
	 * Returns the genome position of the left side of the junction
	 * @return int
	 */
	public int refPosLeft(){
		return this.brl.send;
	}

	/**
	 * Returns the genome position of the right side of the junction.
	 * @return int
	 */
	public int refPosRight(){
		return this.brr.sstart;
	}
	
	/*
	 * 		Read	=============================================>
	 * 		Pos,Pos; no invert
	 * 		Brl		  ------------------>
	 * 		Brr								-------------------->
	 * 		Neg,Neg; invert
	 * 		Brl		  <------------------
	 * 		Brr								<--------------------
	 */
	
	/**
	 * Collect the genome normalized read position for the left side of this junction in read space.
	 * @return int
	 */
	public int queryPosLeft(){
		return this.invert ? this.brr.qstart : this.brl.qend;
	}
	
	/**
	 * Collect the genome normalized read position for the right side of this junction in read space.
	 * @return int
	 */
	public int queryPosRight(){
		return this.invert ? this.brl.qend : this.brr.qstart;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((a == null) ? 0 : a.hashCode());
		result = prime * result + ((b == null) ? 0 : b.hashCode());
		result = prime * result + ((ol == null) ? 0 : ol.hashCode());
		result = prime * result + ((or == null) ? 0 : or.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (!(obj instanceof Junction)) {
			return false;
		}
		Junction other = (Junction) obj;
		if (a == null) {
			if (other.a != null) {
				return false;
			}
		} else if (!a.equals(other.a)) {
			return false;
		}
		if (b == null) {
			if (other.b != null) {
				return false;
			}
		} else if (!b.equals(other.b)) {
			return false;
		}
		if (ol != other.ol) {
			return false;
		}
		if (or != other.or) {
			return false;
		}
		return true;
	}
}
