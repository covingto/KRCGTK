package org.bcm.hgsc.cancer.bed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.utils.RangUtils;

public class BedTools {
	public static final Comparator<BEDRegion> bedsorter = new Comparator<BEDRegion>(){

		@Override
		public int compare(BEDRegion arg0, BEDRegion arg1) {
			int chrcompare = arg0.getSequence().compareTo(arg1.getSequence());
			if (chrcompare == 0){
				return arg0.getStart() - arg1.getStart();
			} else {
				return chrcompare;
			}
		}
		
	};
	
	public static final Comparator<BEDRegion> bedintersected = new Comparator<BEDRegion>(){

		@Override
		public int compare(BEDRegion o1, BEDRegion o2) {
			int chrcompare = o1.getSequence().compareTo(o2.getSequence());
			if (chrcompare == 0){
				if (RangUtils.intersects(o1.getStart(), o1.getEnd(), o2.getStart(), o2.getEnd())){
					return 0;
				} else {
					return o1.getStart() - o2.getStart();
				}
			} else {
				return chrcompare;
			}
		}
		
	};
	
	public static List<BEDRegion> processBedRegions(String chrom, String bedfilename) throws IOException{
		File bedfile = new File(bedfilename);
		if ( ! bedfile.canRead() ){
			throw new IOException("Can't read bed file");
		}
		return processBedRegions(chrom, bedfile);
	}
	
	public static List<BEDRegion> processBedRegions(String chrom,
			File bed) {
		List<BEDRegion> regions = new ArrayList<BEDRegion>();
		// parse a supposed bed file and return a list of BEDRegions.  
		if (chrom == null){
			// no chromosome present, see if there is a bedfile to use
			if (bed == null){
				return null;
			} else {
				if (bed.canRead()){
					try {
						System.out.println("# Reading BED file, no chromosome input.");
						BufferedReader reader = new BufferedReader(new FileReader(bed));
						String line = reader.readLine();
						while (line != null){
							// process the line and add a record to the LinkedList.
							String[] lsplit = StringUtils.split(line, "\t");
							regions.add(new BEDRegion(lsplit[0], Integer.parseInt(lsplit[1]), Integer.parseInt(lsplit[2])));
							line = reader.readLine();
						}
						reader.close();
					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						System.out.println("WARNING; can't read bed file " + bed.getAbsolutePath());
						return null;
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						System.out.println("WARNING; can't read bed file " + bed.getAbsolutePath());
						return null;
					}
				} else {
					System.out.println("WARNING; can't read bed file " + bed.getAbsolutePath());
					return null;
				}
			}
		} else {
			if (bed == null){
				// we have a chr, but no positions in the bed, so we just have to make a single BEDRegion and return that
				regions.add(new BEDRegion(chrom, 0, 0));
			} else {
				if (bed.canRead()){
					try {
						System.out.println("# Reading BED file, " + chrom + ".");
						BufferedReader reader = new BufferedReader(new FileReader(bed));
						String line = reader.readLine();
						System.out.println("Example BED line " + line);
						while (line != null){
							// process the line and add a record to the LinkedList.
							String[] lsplit = StringUtils.split(line, "\t");
							//System.out.println(lsplit[0] + " ? " + chrom + ": " + (lsplit[0].equals(chrom)));
							if (lsplit[0].equals(chrom)){
								regions.add(new BEDRegion(lsplit[0], Integer.parseInt(lsplit[1]), Integer.parseInt(lsplit[2])));
							}
							line = reader.readLine();
						}
						reader.close();
					} catch (FileNotFoundException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						System.out.println("WARNING; can't read bed file " + bed.getAbsolutePath());
						return null;
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						System.out.println("WARNING; can't read bed file " + bed.getAbsolutePath());
						return null;
					}
				} else {
					System.out.println("WARNING; can't read bed file " + bed.getAbsolutePath());
					return null;
				}
			}
		}
		System.out.println("Found " + regions.size() + "BED regions.");
		System.out.println("Sorting bed file");
		Collections.sort(regions, bedsorter);
		return regions;
		
	}
	
	public static int binaryBedHit(List<? extends BEDRegion> regions, BEDRegion query, int imin, int imax){

		if ( imax < imin ){ return -imin; } // return the negative of the location in the list where the variant would be
		int imid = (imin + imax) / 2;
		BEDRegion s = regions.get(imid);

		// test if the key intersects the s
		int hit = bedintersected.compare(s, query);
		if ( hit > 0 ){
			return binaryBedHit(regions, query, imin, imid - 1);
		} else if ( hit < 0 ){
			return binaryBedHit(regions, query, imid + 1, imax);
		} else {
			// key has been found
			return imid;
		}
	}
	
	public static void sort(List<? extends BEDRegion> regions){
		Collections.sort(regions, bedsorter);
	}
}
