package org.bcm.hgsc.utils;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedWriter;
import java.io.IOException;

public class PairedReads {
	final SAMRecord read1;
	final SAMRecord read2;
	
	public PairedReads(SAMRecord samRecord, SAMRecord queryMate) {
		read1 = samRecord;
		read2 = queryMate;
	}
	
	public SAMRecord getRead1(){
		return this.read1;
	}
	
	public SAMRecord getRead2(){
		return this.read2;
	}
	
	public void write(BufferedWriter writer) throws IOException{
		writer.write(read1.getSAMString());
		writer.write(read2.getSAMString());
	}

}
