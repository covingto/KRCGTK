package org.bcm.hgsc.cancer.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class StreamingWalkAcceptor extends WalkAcceptor{
	protected final BufferedWriter writer;
	
	public StreamingWalkAcceptor(File outfile) throws IOException{
		writer = new BufferedWriter(new FileWriter(outfile));
	}
	
	@Override
	public void add(String chr, long posa, long posb, double value) throws IOException{
		writer.write(chr + "\t" + posa + "\t" + posb + "\t" + value + "\n");
	}
	
	public void close() throws IOException{
		this.writer.close();
	}

}
