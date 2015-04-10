package org.bcm.hgsc.cancer.virus;

import htsjdk.samtools.reference.ReferenceSequence;

import java.util.Arrays;

import com.skjegstad.utils.BloomFilter;

public class BloomTrainer implements Runnable {
	
	private final ReferenceSequence ref;
	private final BloomFilter<String> filter;
	
	BloomTrainer(ReferenceSequence ref, BloomFilter<String> filter){
		this.ref = ref;
		this.filter = filter;
	}
	
	@Override
	public void run() {
		System.out.println("Starting " + ref.getName());
		byte[] refbytes = ref.getBases();
		for (int o = 0; o < refbytes.length - 25; o++){
			final byte[] subarray = Arrays.copyOfRange(refbytes, o, o + 25);
			synchronized (filter){
				filter.add(subarray);
			}
		}
		System.out.println("Processed " + ref.getName());
	}

}
