package org.bcm.hgsc.cancer.bed;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class BedGraph {
	private class GraphEntry {
		final String chr;
		final long posa;
		final long posb;
		final double value;
		
		public GraphEntry(String chr, long posa, long posb, double value){
			this.chr = chr;
			this.posa = posa;
			this.posb = posb;
			this.value = value;
		}
	}
	

	private boolean stream = false;
	private String description = "";
	private String visibility = "full";
	private int[] color = { 200, 100, 0 };
	private int[] altColor = { 0, 100, 200 };
	private String name = "";
	private List<GraphEntry> entries = new ArrayList<GraphEntry>();
	private Set<String> chrs = new HashSet<String>();
	private BufferedWriter writer = null;
	
	public BedGraph(){
		this("", "");
	}

	public BedGraph(String name, String description) {
		this.name = name;
		this.description = description;
	}
	
	public double stdev(){ 
		double avg = this.average();
		double sqmeandiff = 0;
		for (int i = 0; i < entries.size(); i++){
			sqmeandiff += Math.pow((entries.get(i).value - avg), 2);
		}
		return Math.sqrt(sqmeandiff / entries.size());
	}
	
	// averages
	public Double average(){
		double total = 0;
		int count = 0;
		for (GraphEntry ge : entries){
			total += ge.value;
			count += 1;
		}
		return total / count;
	}
	
	public Map<String, Double> chrAverage(){
		Map<String, Double> totals = new HashMap<String, Double>();
		Map<String, Integer> counts = new HashMap<String, Integer>();
		Map<String, Double> averages = new HashMap<String, Double>();
		
		for (GraphEntry ge : entries){
			if (!totals.containsKey(ge.chr)){
				totals.put(ge.chr, 0.0);
			}
			if (!counts.containsKey(ge.chr)){
				counts.put(ge.chr, 0);
			}
			totals.put(ge.chr, totals.get(ge.chr) + ge.value);
			counts.put(ge.chr, counts.get(ge.chr) + 1);
		}
		
		for (String k : totals.keySet()){
			averages.put(k, totals.get(k) / counts.get(k));
		}
		return averages;
	}
	
	// adding samples
	public void add(String chr, long posa, long posb, double value) throws IOException{
		if (this.stream && this.writer != null){
			writer.write(chr + "\t" + posa + "\t" + posb + "\t" + value);
			writer.newLine();
		} else {
			chrs.add(chr);
			entries.add(new GraphEntry(chr, posa, posb, value));
		}
	}
	
	public List<GraphEntry> getEntries(){ 
		return this.entries;
	}
	
	// writers
	public void write(File filename) throws IOException{
		FileWriter fw = new FileWriter(filename);
		BufferedWriter writer = new BufferedWriter(fw);
		write(writer);
		writer.close();
		fw.close();
	}
	
	public void write(BufferedWriter writer) throws IOException{
		writer.write("# average=" + this.average() + " stdev=" + this.stdev());
		writer.newLine();
		writer.write("track type=bedGraph name=\"" + this.name + "\" description=\"" + this.description + 
				"\" color=" + this.color[0] + "," + this.color[1] + "," + this.color[2] + 
				" altColor=" + this.altColor[0] + "," + this.altColor[1] + "," + this.altColor[2] + 
				" visibility=" + this.visibility);
		for (GraphEntry ge : entries){
			writer.newLine();
			writer.write(ge.chr + "\t" + ge.posa + "\t" + ge.posb + "\t" + ge.value);
		}
		writer.newLine();
	}
	
	// modifiers
	public void setVisibility(String v){
		this.visibility = v;
	}

	public void setColor(int r, int g, int b){
		this.color[0] = r;
		this.color[1] = g;
		this.color[2] = b;
	}
	
	public void setAltColor(int r, int g, int b){
		this.altColor[0] = r;
		this.altColor[1] = g;
		this.altColor[2] = b;
	}

	
	public void setStreaming(File mmFile) throws IOException {
		this.stream = true;
		this.writer = new BufferedWriter(new FileWriter(mmFile));
	}
	
	public void close() throws IOException{
		if (this.writer != null){
			this.writer.close();
		}
	}
}
