package org.bcm.hgsc.cancer.utils;

import java.io.File;
import java.io.IOException;

public class LimitingStreamingWalkAcceptor extends StreamingWalkAcceptor {
	private final double minval;
	private final double maxval;
	
	public LimitingStreamingWalkAcceptor(File outfile, double minval, double maxval) throws IOException {
		super(outfile);
		this.minval = minval;
		if (maxval > 0) { this.maxval = maxval; }
		else { this.maxval = Double.MAX_VALUE; }
	}
	
	@Override
	public void add(String chr, long posa, long posb, double value) throws IOException{
		if (value >= this.minval && value <= this.maxval){
			writer.write(chr + "\t" + posa + "\t" + posb + "\t" + value + "\n");
		}
	}
}
