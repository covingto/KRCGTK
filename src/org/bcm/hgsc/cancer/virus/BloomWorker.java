package org.bcm.hgsc.cancer.virus;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

import java.io.File;

import org.bcm.hgsc.utils.BAMInterface;
import org.bcm.hgsc.utils.Utils;

public class BloomWorker extends BloomQuery implements Runnable  {
	
	private final SynchronizedSamWriter writer;
	private final SamReader bam;
	
	BloomWorker(String filterFile, String bamFile, SynchronizedSamWriter w) throws Exception{
		// here we extract the file
		super(filterFile);
		writer = w;
		bam = new BAMInterface(new File(bamFile), null, null).getSamfilereader();
	}
	
	@Override
	public void run() {
		SAMRecord rec;
		try {
			while( (rec = ViralBloomFilter.samRecords.take()) != null){
				String seq = rec.getReadString();
				//byte[] bytes = rec.getReadBases();
				String revcomp = Utils.reverseComplement(seq);
				if (test(seq) || test(revcomp)){
					// get the mate
					ViralBloomFilter.hitCount++;
					SAMRecord mate = bam.queryMate(rec);
					synchronized (writer){
						writer.writeHit(rec);
						ViralBloomFilter.matesToGet.add(mate);
						//writer.writeMate(mate);
					}
				}
			}
			ViralBloomFilter.samRecords.add(null); // propogate the signal
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	

}
