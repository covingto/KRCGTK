package org.bcm.hgsc.cancer.virus;

import java.io.IOException;
import java.util.Arrays;

import com.skjegstad.utils.BloomFilter;

public class BloomQuery {
	protected final BloomFilter<String> filter;
	
	BloomQuery(String filterFile) throws IOException, ClassNotFoundException{
		filter = ViralBloomFilter._deserialize(filterFile);
	}
	
	/**
	 * Tests if a sequence might be in the set.  This is done by splitting the sequence into three 25mer sections
	 * and testing each.  If 2 of three return true then the test returns true.
	 * @param s
	 */
	public boolean test(String s){
		byte[] bytes = s.getBytes();
		int len = bytes.length;
		int offset = (len - (len % 25)) / 25;
		int m = 0;
		int o1 = offset;
		int o2 = (offset * 2) + 25;
		int o3 = (offset * 3) + (25 * 2);
		byte[] t1 = Arrays.copyOfRange(bytes, o1, o1 + 25);
		byte[] t2 = Arrays.copyOfRange(bytes, o2, o2 + 25);
		byte[] t3 = Arrays.copyOfRange(bytes, o3, o3 + 25);
		
		if (filter.contains(t1)){ m++; }
		if (filter.contains(t2)){ m++; }
		if (m == 2){ return true; } // this is just a good time to do this since it saves us one check
		if (m == 0){ return false; } // can never make it
		if (filter.contains(t3)){ m++; }
		return m >= 2;
	}
}
