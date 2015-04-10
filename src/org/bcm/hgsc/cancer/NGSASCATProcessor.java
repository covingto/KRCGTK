package org.bcm.hgsc.cancer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;

public class NGSASCATProcessor {
	private final File pairsfile;
	private final File ascatdir;
	private final File mafdir;
	private final File segdir;
	private String[][] pairsdata;
	
	NGSASCATProcessor(File pairsfile, File ascatdir, File mafdir, File segdir) throws IOException{
		this.pairsfile = pairsfile;
		this.ascatdir = ascatdir;
		this.mafdir = mafdir;
		this.segdir = segdir;
		this.pairsdata = CancerUtils.readPairsFile(this.pairsfile);
	}
	
//	NGSASCATProcessor(File projdir){
//		// find or create all required directories from a standard project direcory
//		FilenameFilter varscanFilter = new CancerUtils.DirectoryFilter("varscan");
//		FilenameFilter copynumberFilter = new CancerUtils.DirectoryFilter("copynumber");
//		
//		
//		File varscandir = projdir.listFiles(varscanFilter)[0];
//		File copynumberdir = varscandir.listFiles(copynumberFilter)[0];
//		this.segdir = copynumberdir;
//	}
	
	public void readData(){
		// read data from a series of maf and seg files to build the data structure.
		HashMap<String, SegmentationData> segdata = new HashMap<String, SegmentationData>();
		HashMap<String, MAFData> mafdata = new HashMap<String, MAFData>(); 
		
		class MafFinder implements FilenameFilter{
			@Override
			public boolean accept(File dir, String name) {
				if (name.endsWith("mafplus")) return true;
				return false;
			}
		}
		
		class SegFinder implements FilenameFilter{
			@Override
			public boolean accept(File dir, String name) {
				if (name.endsWith("seg")) return true;
				return false;
			}
		}
		
		File[] segfiles = segdir.listFiles(new SegFinder());
		File[] maffiles = mafdir.listFiles(new MafFinder());
		
		for (final File segfile : segfiles){
			BufferedReader segreader = null;
			String line;
			try {
				segreader = new BufferedReader(new FileReader(segfile));
				SegmentationData segments = new SegmentationData();
				while ((line = segreader.readLine()) != null){
					if (line.contains("ID")) continue; // this is the header line
					String[] lsplit = line.split("\t");
					segments.add(lsplit[1], Integer.parseInt(lsplit[2]), Integer.parseInt(lsplit[3]), Integer.parseInt(lsplit[4]), Float.parseFloat(lsplit[5]));
				}
				segreader.close();
				
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} finally {
				if (segreader != null){
					try {
						segreader.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		
		
		
	}

}
