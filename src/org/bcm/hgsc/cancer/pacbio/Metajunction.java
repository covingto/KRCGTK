package org.bcm.hgsc.cancer.pacbio;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.cancer.sv.Junction;
import org.bcm.hgsc.cancer.utils.Chromosome;

public class Metajunction {
	private List<Junction> junctions = new ArrayList<Junction>();
	private List<Integer> posl = new ArrayList<Integer>();
	private List<Integer> posr = new ArrayList<Integer>();
	private List<String> reads = new ArrayList<String>();
	private List<Integer> qposl = new ArrayList<Integer>();
	private List<Integer> qposr = new ArrayList<Integer>();
	private final Chromosome chrl;
	private final Chromosome chrr;
	private final int buffer;
	
	public Metajunction(Junction seed, int buffer){
		this.chrl = seed.a.chr;
		this.chrr = seed.b.chr;
		this.buffer = buffer;
		this.add(seed);
	}

	private static long sum(List<Integer> list) {
	     long sum= 0; 
	     for (Integer i:list)
	         sum = sum + i;
	     return sum;
	}
	
	public boolean fits(Junction test) {
		Chromosome tchrl = test.a.chr;
		Chromosome tchrr = test.b.chr;
		Integer tposl = test.a.getStart();
		Integer tposr = test.b.getStart();
		
		if ( ! chrl.equals(tchrl) || ! chrr.equals(tchrr) ){ return false; }
		long meanl = Math.round(sum(posl) / (double) posl.size());
		long meanr = Math.round(sum(posr) / (double) posr.size());
		if ( Math.abs(tposl - meanl) < buffer && Math.abs(tposr - meanr) < buffer ){
			return true;
		}
		return false;
	}

	public void add(Junction test) {
		this.junctions.add(test);
		this.posl.add(test.a.getStart());
		this.posr.add(test.b.getStart()); 
		this.reads.add(test.a.read + "(" + test.brl.qstart + "-" + test.brl.qend + ">" + test.brr.qstart + "-" + test.brr.qend + ";" + test.invert + ")");
		this.qposl.add(test.queryPosLeft());
		this.qposr.add(test.queryPosRight());
	}

	public int size() {
		return this.junctions.size();
	}
	
	public String ori(){
		Set<String> seen = new HashSet<String>();
		for ( Junction j : junctions ){
			seen.add(j.ol + "," + j.or);
		}
		return StringUtils.join(seen, "|");
	}
	
	public String leftStatus(){
		Set<String> status = new HashSet<String>();
		for ( Junction j : junctions ){
			status.add(j.brl.unique.toString());
		}
		return StringUtils.join(status, ",");
	}
	
	public String rightStatus(){
		Set<String> status = new HashSet<String>();
		for ( Junction j : junctions ){
			status.add(j.brr.unique.toString());
		}
		return StringUtils.join(status, ",");
	}
	
	public long gapDist(){
		long dist = 0;
		for (Junction j : this.junctions){
			dist += j.readGap();
		}
		return Math.round(dist / (double) this.junctions.size());
	}
	
	@Override
    public String toString(){
		// sanity check that we are writing positive junctions
		// at some point it was noticed that negative junctions were reported, at the approximate location of the proper junction (just in negative space)
		// may be commented later with a report
		
		// calculate the mapping location, these are the locations of the junctions in genome normalized space
		long meanl = Math.round(sum(this.posl) / (double) this.posl.size());
		long meanr = Math.round(sum(this.posr) / (double) this.posr.size());
		long meanreadgap = gapDist();
		String junctiontype = "UNK";
		String meanrefgapstr = "NA";
		if (this.chrl.equals(this.chrr)){
			// this is an interchromosomal event, can now be INS, DEL, ITX
			
			// check to see what the orientations are, if they are equal (POS, POS) then this can be INS or DEL
			boolean itx = false;
			for ( Junction j : junctions ){
				if (!j.ol.equals(j.or)){
					itx = true;
					break;
				}
			}
			if (itx){
				junctiontype = "ITX";
				meanrefgapstr = "NA";
			} else {
				long meanrefgap = meanr - meanl;
				meanrefgapstr = String.valueOf(meanrefgap);
				if ((meanreadgap <= 0 && meanrefgap <= 0) || (meanreadgap > 0 && meanrefgap > 0)){
					if (meanreadgap < meanrefgap){
						/*
						 *   read  -----------------------------------------------------
						 *   map1  =====================================
						 *   map2                    ===================================
						 *   ref   -----------------------------------------------------
						 *   map1  ===========================
						 *   map2                            ===========================
						 *                          Deletion in the read, but complex
						 */
						junctiontype = "CDEL"; // complex deletion
					} else {
						junctiontype = "CINS"; // complex insertion
					}
				} else if (meanreadgap <= 0 && meanrefgap > 0){
					junctiontype = "DEL"; // these are always deletions
				} else if (meanreadgap > 0 && meanrefgap <= 0){
					junctiontype = "INS"; // these are always insertions
				}
			}
		} else {
			junctiontype = "CTX";
		}
		// see if the blast rows are ambiguous
		String lstatus = leftStatus();
		String rstatus = rightStatus();
		
		if (meanl < 0 | meanr < 0){
			System.err.println("Negative position for metajunction.");
			System.err.println("Posl; " + meanl + " Posl size; " + posl.size() + " Posr; " + meanr + " Posr size; " + posr.size());
			System.err.println(StringUtils.join(posl, ","));
			System.err.println(StringUtils.join(posr, ","));
		}
		
		return this.junctions.size() + "\t" + 
				// print the genomic positions
				this.chrl + "\t" + 
				meanl + "\t" + 
				this.chrr + "\t" + 
				meanr +  "\t" +
				// print the junction type
				junctiontype + "\t" + 
				// print the gap information for this junction
				meanrefgapstr + "\t" + 
				meanreadgap + "\t" +
				// print the mapping ori
				this.ori() + "\t" +
				// print the read data
				StringUtils.join(reads, ",") + "\t" +
			    // print some information about the reads themselves
				lstatus + "\t" +
				rstatus
				;
	}
}
