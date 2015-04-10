package org.bcm.hgsc.cancer.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.cancer.MAFData;

public class Mafstruct {
	public static Comparator<MAFData> varcomparator = new Comparator<MAFData>(){

		@Override
		public int compare(MAFData arg0, MAFData arg1) {
			if ( arg0.get("Chromosome") == arg1.get("Chromosome") ){
				if ( Integer.parseInt(arg0.get("Start_position")) == Integer.parseInt(arg1.get("Start_position"))){
					if ( Integer.parseInt(arg0.get("End_position")) == Integer.parseInt(arg1.get("End_position")) ){
						if ( arg0.get("Tumor_Seq_Allele2").equals(arg1.get("Tumor_Seq_Allele2")) ){
							if ( arg0.get("Tumor_Sample_Barcode").equals(arg1.get("Tumor_Sample_Barcode")) ){
								if ( arg0.get("Algorithm").equals(arg1.get("Algorithm")) ){
									return 0;
								} else {
									return arg0.get("Algorithm").compareTo(arg1.get("Algorithm"));
								}
							} else {
								return arg0.get("Tumor_Sample_Barcode").compareTo(arg1.get("Tumor_Sample_Barcode"));
							}
						} else {
							return arg0.get("Tumor_Seq_Allele2").compareTo(arg1.get("Tumor_Seq_Allele2"));
						}
					} else {
						return Integer.parseInt(arg0.get("End_position")) - Integer.parseInt(arg1.get("End_position"));
					}
				} else {
					return Integer.parseInt(arg0.get("Start_position")) - Integer.parseInt(arg1.get("Start_position"));
				}
			} else {
				return arg0.get("Chromosome").compareTo(arg1.get("Chromosome"));
			}
		}
		
	};
	
	private List<MAFData> variants = new ArrayList<MAFData>();
	private Set<String> subjects = new HashSet<String>();
	private List<String> header = null;
	
	
	public Mafstruct(File f) throws IOException{
		// read in data from a file and store it in the structure
		
		FileReader fr = new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		header = Arrays.asList(br.readLine().split("\t"));
		// note we need to handle the failed reason this could be in the maf header or not
		
		String line;
		List<String> lsplit;
		System.out.println("Reading variants");
		int varcount = 0;
		while ((line = br.readLine()) != null) {
		   // process the line.
			varcount++;
			lsplit = Arrays.asList(line.split("\t"));
			MAFData newmd = new MAFData(header, lsplit);
			variants.add(newmd);
			subjects.add(newmd.get("Tumor_Sample_Barcode"));
			if (varcount%1000 == 0){
				System.out.println("Added " + varcount + "variants");
			}
		}
		br.close();
		Collections.sort(variants, varcomparator);
	}
	
	public void sort(){
		sort(varcomparator);
	}
	
	public void sort(Comparator<MAFData> comp){
		Collections.sort(variants, comp);
	}

	public int length() {
		return variants.size();
	}
	
	public MAFData get(int i) {
		return variants.get(i);
	}
	
	public List<MAFData> getPassing() {
		List<MAFData> result = new ArrayList<MAFData>();
		
		for (MAFData md : variants) {
			if (md.get("Failed_Reason").equals("Pass")){
				result.add(md);
			}
		}
		
		return result;
	}
	
	public List<MAFData> getFailing(){
		List<MAFData> result = new ArrayList<MAFData>();
		
		for (MAFData md : variants) {
			if ( ! md.get("Failed_Reason").equals("Pass")){
				result.add(md);
			}
		}
		
		return result;
	}

	
	private void writeVariants(List<MAFData> md, File outfile){
		BufferedWriter writer = null;
		List<String> outheader = new ArrayList<String>(header);
		if ( ! outheader.contains("Failed_Reason") ) {
			outheader.add("Failed_Reason");
		}
		try {
			writer = new BufferedWriter( 
					new OutputStreamWriter(
							new FileOutputStream(outfile), "utf-8"
							)
					);
		
			writer.write(StringUtils.join(outheader, "\t"));
			writer.newLine();
			
			for ( MAFData v : md ) {
				writer.write(v.toString(outheader));
				writer.newLine();
			}
		} catch (IOException ex) {
			ex.printStackTrace();
		} finally{
			try {writer.close(); } catch (Exception ex) {}
		}
	}
	

	public void writePassing(String somaticPassOut) {
		writePassing(new File(somaticPassOut));
	}
	
	public void writePassing(File somaticPassOut) {
		writeVariants(getPassing(), somaticPassOut);
	}
	
	public void writeFailing(String somaticFailOut) {
		writeFailing(new File(somaticFailOut));
	}
	
	public void writeFailing(File somaticFailOut){
		writeVariants(getFailing(), somaticFailOut);
	}
	
	public List<MAFData> getVariants(){
		return variants;
	}
	
	public Set<String> subjects() {
		return subjects;
	}
	
	/**
	 *   Returns a List of MAFData pertaining to the specified subject.  If the subject isn't in the MAF an empty list is returned.
	 */
	public List<MAFData> subsetSubject(String s){
		List<MAFData> result = new ArrayList<MAFData>();
		
		for (MAFData md : variants) {
			if (md.get("Tumor_Sample_Barcode").equals(s)){
				result.add(md);
			}
		}
		return result;
	}
	
}
