package org.bcm.hgsc.cancer;

import java.util.LinkedList;

/**
 * SegmentationData class
 * 
 * <p>
 * This class represents a data structure which contains rows of segmentation data.  The class provides accessory functions
 * to interact with this data.
 * </p>
 * @author covingto
 *
 */
public class SegmentationData {
	class Segment {
		String chrom;
		int start;
		int end;
		int nummarks;
		float segmean;
		
		Segment(String chrom, int start, int end, int nummarks, float segmean){
			this.chrom = chrom;
			this.start = start;
			this.end = end;
			this.nummarks = nummarks;
			this.segmean = segmean;
		}
	}
	
	private LinkedList<Segment> segments = new LinkedList<Segment>();
	
	public void add(String chrom, int start, int end, int nummarks, float segmean){
		segments.add(new Segment(chrom, start, end, nummarks, segmean));
	}
	
	public LinkedList<Segment> getSegments() {
		return segments;
	}
}
