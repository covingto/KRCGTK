package org.bcm.hgsc.cancer.pacbio;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import org.bcm.hgsc.cancer.utils.Chromosome;
import org.bcm.hgsc.utils.Settings;

public class BlastStruct {
	
	private Map<String, BlastGroup> mappings = new HashMap<String, BlastGroup>();
	
	public BlastStruct(File blastres, IndexedFastaSequenceFile reffasta) throws IOException{
		this(blastres, false, false, reffasta);
	}
	
	public BlastStruct(File blastres) throws IOException{
		this(blastres, false, false, null);
	}
	
	public BlastStruct(File blastres, boolean allowrandom, boolean allowun, IndexedFastaSequenceFile reffasta) throws IOException{
		FileReader fr = new FileReader(blastres);
		BufferedReader reader = new BufferedReader(fr);
		String line;
		String[] lsplit;
		int l = 0;
		System.out.println("Beginnign File read");
		System.out.println(Settings.dateFormat.format(new Date()));
		while((line = reader.readLine()) != null) {
			// see if we need to make a new collector
			l++;
			if (Settings.debug && l % 1000 == 0){
				System.out.print("Processed " + l + " lines.\r");
				if (l % 100000 == 0){ System.out.println("Processed " + l + " lines."); }
			}
			lsplit = line.split("\t");
			
			// do some basic read filtering
			if (!allowrandom & lsplit[1].contains("random")) { 
				if (Settings.debug) { System.out.println("Ignoring; " + line); }
				continue; } // we don't allow random mappings
			if (!allowun & lsplit[1].contains("chrUn")) { 
				if (Settings.debug) { System.out.println("Ignoring; " + line); }
				continue; }
			// this is where we generate the blast row object to take advantage of some features
			BlastRow br = new BlastRow(lsplit, reffasta);
			Chromosome chr = new Chromosome(br.schr);
			if (chr.getChrnum() == 0){ 
				if (Settings.debug) { System.out.println("Ignoring; " + line); }
				continue; } // we don't work with non-standard chromosomes.
			// add the group if it doesn't exist
			if (!mappings.containsKey(lsplit[0])){ mappings.put(lsplit[0], new BlastGroup(lsplit[0], Integer.parseInt(lsplit[12]))); }
			
			mappings.get(lsplit[0]).add(br);
		}
		reader.close();
		System.out.println("Reading of file finished");
		System.out.println(Settings.dateFormat.format(new Date()));
		System.out.println("Sorting groups by qstart position, and setting unique status");
		
		// do some group filtering to label reads as redundant or otherwise questionable
		/**
		 * 	Filter based on mapping characteristics to indicate the status with respect to the global alignment pattern.
		 *   <pre>
		 *   READ -----------------------------------------------------------------
		 *   MAP1    -----------------------										UNIQUE
		 *   MAP2                                       -------------------			AMBIGUOUS
		 *   MAP3                                        -------------------		AMBIGUOUS
		 *   MAP4                                              ---					REDUNDANT
		 *   </pre>
		 */
		
		for (BlastGroup bg : mappings.values()){
			bg.sort();
			for (int i = 0; i < (bg.numMappings() - 1); i++){
				// collect mapping i
				BlastRow m1 = bg.getMapping(i);
				for (int j = i + 1; j < (bg.numMappings()); j++){
					// collect mapping j (i + n)
					BlastRow m2 = bg.getMapping(j);
					// check if we excede the end size, these can't be redundant so break to avoid excess computation
					if (m1.qend < m2.qstart){ break; } 
					// calculte m1 overlap pct
					Float m1ol = BlastRow.overlapFrc(m1.qstart, m1.qend, m2.qstart, m2.qend);
					Float m2ol = BlastRow.overlapFrc(m2.qstart, m2.qend, m1.qstart, m1.qend);
					if ( m1ol > 0.90 & m2ol > 0.90 ) { 
						// we don't reassign from redundant as that is the lowest quality of mapping available
						if (m1.unique != BlastRow.REDUNDANT) { m1.unique = BlastRow.AMBIGUOUS; }
						if (m2.unique != BlastRow.REDUNDANT) { m2.unique = BlastRow.AMBIGUOUS; }
					} 
					else if ( m1ol > 0.90 ){ m1.unique = BlastRow.REDUNDANT; }
					else if ( m2ol > 0.90 ){ m2.unique = BlastRow.REDUNDANT; }
					if ( Settings.debug ){
						System.out.println("Compare; " + m1.qstart + " " + m1.qend + " " + 
								m2.qstart + " " + m2.qend + " " + m1ol + " " + m2ol + " " + 
								m1.unique + " " + m2.unique);
					}
				}
			}
		}
		
	}
	
	public Map<String, BlastGroup> getMappings(){
		return mappings;
	}
	
	public int numReads(){
		return this.mappings.size();
	}
	
	public int numAlignments(){
		int n = 0;
		for (String k : this.mappings.keySet()){
			BlastGroup bg = this.mappings.get(k);
			n += bg.mappingLength();
		}
		
		return n;
	}
	
	public void printSummary(){
		System.out.println("Numer of reads: " + this.numReads());
		System.out.println("Number of alignments: " + this.numAlignments());
	}
	
	
}
