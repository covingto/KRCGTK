package org.bcm.hgsc.cancer.virus;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

import java.io.File;

public class SynchronizedSamWriter {
	private final SAMFileWriter writer;
	SynchronizedSamWriter(String outputFile, SamReader inputSam){
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		factory.setMaxRecordsInRam(100);
		SAMFileHeader newHeader = inputSam.getFileHeader().clone();
		newHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		writer = factory.makeBAMWriter(newHeader, false, new File(outputFile));
	}
	
	public void writeHit(SAMRecord rec){
		writer.addAlignment(rec);
	}
	
	public void close(){
		writer.close();
	}
	
	public void writeMate(SAMRecord rec){
		// writing a mate checks that the mate has not already been written as a hit
	}
}
