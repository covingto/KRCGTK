package org.bcm.hgsc.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.bcm.hgsc.cancer.bed.BEDRegion;

public class PairSet {
	//private List<SAMRecord> reads = new ArrayList<SAMRecord>();
	private List<PairedReads> pairs = new ArrayList<PairedReads>();
	private List<SAMRecord> unpaired = new ArrayList<SAMRecord>();
	private int topair;
	private final BEDRegion r;
	
	public PairSet(String bamfile, BEDRegion r) throws IOException{
		this.r = r;
		List<SAMRecord> slicereads = querySet(bamfile, r.getSequence(), r.getStart(), r.getEnd());
		List<SAMRecord> chunkreads = querySet(bamfile, r.getSequence(), r.getStart() - FindPairs.chunkbuffer, r.getEnd() + FindPairs.chunkbuffer);
		
		boolean[] flagged = new boolean[slicereads.size()];
		Arrays.fill(flagged, false);
		
		// fill from the initial slice
		for (int i = 0; i < slicereads.size() - 1; i++){
			if (flagged[i]) continue; // flagged means that we have already put the read into a pair and don't need to do that again
			final String readName = slicereads.get(i).getReadName();
			for (int j = i + 1; j < slicereads.size(); j++){
				if (flagged[j]) continue;
				if (slicereads.get(j).getReadName().equals(readName)){
					this.pairs.add(new PairedReads(slicereads.get(i), slicereads.get(j)));
					flagged[i] = true;
					flagged[j] = true;
					break; // found the mate, stop looking for more
				}
			}
			if (!flagged[i]){ // look up in the chunk
				boolean firstofpair = slicereads.get(i).getFirstOfPairFlag();
				for (int j = 0; j < chunkreads.size(); j++){
					SAMRecord next = chunkreads.get(j);
					if (next.getReadName().equals(readName) && next.getFirstOfPairFlag() != firstofpair){
						this.pairs.add(new PairedReads(slicereads.get(i), chunkreads.get(j)));
						flagged[i] = true;
						break;
					}
				}
			}
			if (!flagged[i]){ // need to add this to the unpaired list
				this.unpaired.add(slicereads.get(i));
			}
		}
		this.topair = this.unpaired.size(); // we will decrease this as we make pairs
	}
	
	public boolean pairHit(SAMRecord e){
		// if everything is already paired, we don't need to pair anything just return false to let the other guys have a shot at it
		if (topair < 1) return false;
		for (int i = 0; i < unpaired.size(); i++){
			final SAMRecord u = unpaired.get(i);
			if (u.getReadName().equals(e.getReadName()) && u.getFirstOfPairFlag() != e.getFirstOfPairFlag()){ // yay we found a pair
				topair = topair - 1;
				pairs.add(new PairedReads(u, e));
				return true;
			}
		}
		// if we don't find anything then we didn't hit a pair
		return false;
	}
	
	public void writePairs() throws IOException{
		BufferedWriter writer = new BufferedWriter(new FileWriter("reads" + this.r.getSequence() + "_" + this.r.getStart() + "_" + this.r.getEnd() + ".fp.sam")); 
		for (PairedReads p : pairs){
			p.write(writer);
		}
		writer.close();
	}
	

	private List<SAMRecord> querySet(String bamfile, String chr, int start, int end) throws IOException{
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		SamReader inputSam = SamReaderFactory.makeDefault().open(new File(bamfile));
		SAMRecordIterator sri = inputSam.query(chr, start, end, false);
		// get the reads for the slice
		while (sri.hasNext()){
			final SAMRecord sr = sri.next();
			// for these purposes we only are interested in paired reads
			if (!sr.getReadPairedFlag()) continue;
			reads.add(sr);
		}
		sri.close();
		inputSam.close();
		return reads;
	}
}
