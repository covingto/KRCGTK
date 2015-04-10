package org.bcm.hgsc.cancer.pacbio;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import org.bcm.hgsc.cancer.bed.BEDRegion;
import org.bcm.hgsc.cancer.bed.BedTools;
import org.bcm.hgsc.cancer.utils.Orientation;

public class BlastGroup {
	public static Comparator<BlastRow> qstartComp = new Comparator<BlastRow>() {

		@Override
		public int compare(BlastRow arg0, BlastRow arg1) {
			// TODO Auto-generated method stub
			return arg0.qstart.compareTo(arg1.qstart);
		}
		
	};
	
	private List<BlastRow> mappings = new ArrayList<BlastRow>();
	private final Integer length;
	
	public BlastGroup(String q, Integer len){
		this.length = len;
	}
	
	public void add(BlastRow r){
		mappings.add(r);
	}
	
	public Integer length(){
		return length;
	}
	
	public Integer numMappings(){
		return mappings.size();
	}
	
	public Integer numPassingMappings(){
		Integer i = 0;
		for (BlastRow br : mappings){
			if (br.unique == BlastRow.AMBIGUOUS | br.unique == BlastRow.UNIQUE){
				i++;
			}
		}
		return i;
	}
	
	public void sort(){
		Collections.sort(mappings, qstartComp);
	}
	
	public BlastRow getMapping(Integer i){
		return mappings.get(i);
	}
	
	public List<BlastRow> getUnique(){
		List<BlastRow> result = new ArrayList<BlastRow>();
		for (BlastRow br : mappings){
			if (br.unique == BlastRow.UNIQUE) {
				result.add(br);
			}
		}
		return result;
	}

	public void write(BufferedWriter writer) throws IOException {
		// TODO Auto-generated method stub
		for (BlastRow br : mappings){
			writer.write(br.toString());
			writer.newLine();
		}
	}
	
	public static void filterUnmapped(BlastStruct reads, File unmapped) throws IOException{
		// move unmapped reads to the unmapped file
		FileWriter fw = new FileWriter(unmapped);
		BufferedWriter writer = new BufferedWriter(fw);
		Map<String, BlastGroup> mappings = reads.getMappings();
		List<String> toRemove = new ArrayList<String>();
		for (String k : mappings.keySet()){
			if (isBlastUnmapped(mappings.get(k))){
				mappings.get(k).write(writer);
				toRemove.add(k);
			}
		}
		writer.close();
		
		System.out.println(toRemove.size() + " reads were excluded from break analysis.");
		System.out.println((mappings.size() - toRemove.size()) + " reads continue for break analysis.");
		
		// now that itteration is over we actually remove the data
		for (String k : toRemove){
			mappings.remove(k);
		}
	}
	
	public static void filterFullyMapped(BlastStruct reads, File fullymapped) throws IOException{
		// move unmapped reads to the unmapped file
		FileWriter fw = new FileWriter(fullymapped);
		BufferedWriter writer = new BufferedWriter(fw);
		Map<String, BlastGroup> mappings = reads.getMappings();
		List<String> toRemove = new ArrayList<String>();
		for (String k : mappings.keySet()){
			if (isBlastFullyMapped(mappings.get(k), 0.9)){
				mappings.get(k).write(writer);
				toRemove.add(k);
			}
		}
		writer.close();

		System.out.println(toRemove.size() + " reads were fully mapped.");
		System.out.println((mappings.size() - toRemove.size()) + " reads continue analysis.");

		// now that itteration is over we actually remove the data
		for (String k : toRemove){
			mappings.remove(k);
		}
	}
	
	private static boolean isBlastFullyMapped(BlastGroup g, double d){
		for (int i = 0; i < g.mappings.size(); i++){
			BlastRow br = g.mappings.get(i);
			if ( ((double) Math.abs(br.qend - br.qstart) / (double) g.length()) > d){ return true; }
		}
		return false;
	}
	
	private static boolean isBlastUnmapped(BlastGroup g) {
		// TODO Auto-generated method stub
		// Integer readlen = g.length();
		if (g.numPassingMappings() < 2){
			// there is only one mapping for this group, so we just let it go
			// System.out.println("Number of mappings less than 2, unmapped");
			return true;
		}
		if (g.getUnique().size() < 1){
			// no unique values, the entire thing is ambiguous
			// System.out.println("Number of unique is less than 1, unmapped");
			return true;
		}
		return false;
	}
	
	public int mappingLength(){
		return this.mappings.size();
	}

	public static void filterBed(BlastStruct reads, File bedfile){
		filterBed(reads, bedfile, 0);
	}
	
	public static void filterBed(BlastStruct reads, File bedfile, int buffer) {
		List<BEDRegion> regions = BedTools.processBedRegions(null, bedfile);
		filterBed(reads, regions, buffer);
	}

	private static boolean filterBedGroup(BlastGroup blastGroup,
			List<? extends BEDRegion> regions, int buffer) {
		for (BlastRow br : blastGroup.mappings ){
			// is the row in the region?
			BEDRegion thisBed;
			if (br.o == Orientation.POS){
				thisBed = new BEDRegion(br.schr, br.sstart - buffer, br.send + buffer);
			} else {
				thisBed = new BEDRegion(br.schr, br.send - buffer, br.sstart + buffer);
			}
			if (BedTools.binaryBedHit(regions, thisBed, 0, regions.size()) >= 0){
				return false;
			}
		}
		return true;
	}
	

	public static void filterBed(BlastStruct reads, List<? extends BEDRegion> regions,
			int buffer) {
		List<String> toRemove = new ArrayList<String>();
		for (String k : reads.getMappings().keySet() ){
			if (filterBedGroup(reads.getMappings().get(k), regions, buffer)){
				toRemove.add(k);
			}
		}
		for (String k : toRemove){
			reads.getMappings().remove(k);
		}
	}
}
