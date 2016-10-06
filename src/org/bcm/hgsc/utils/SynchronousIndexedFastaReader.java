package org.bcm.hgsc.utils;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class SynchronousIndexedFastaReader {
	
	public static synchronized ReferenceSequence getSubsequenceAt(IndexedFastaSequenceFile fasta, String contig, long start, long stop ) {
		int attempts = 0;
		while (attempts < 10){
			try {
				// fasta.reset(); // let's reset twice just to make sure that we got it :)
				ReferenceSequence result = fasta.getSubsequenceAt(contig, start, stop);
				// fasta.reset();
				return result;
			} catch(SAMException ex){
				attempts += 1;
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					continue;
				}
			}
		}
		throw new SAMException("Unable to load " + contig + "(" + start + ", " + stop + ")");
	}
	
	public static byte[] getBytesAt(IndexedFastaSequenceFile fasta, String contig, long start, long stop ){
		ReferenceSequence seq = getSubsequenceAt(fasta, contig, start, stop);
		return seq.getBases();
	}
	
	public static String getSeqStringAt(IndexedFastaSequenceFile fasta, String contig, long start, long stop ){
		return new String(getBytesAt(fasta, contig, start, stop));
	}
}
