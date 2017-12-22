package org.bcm.hgsc.utils;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

public class BaseQualityAdjust {

	public static void main(String[] args) throws ParseException, IOException {
		// TODO Auto-generated method stub
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("i", true, "input bam file name");
		options.addOption("o", true, "output bam file name");
		options.addOption("q", true, "original quality of the bam");
		options.addOption("n", true, "new quality of the bam");
		CommandLine line = parser.parse(options, args);
		File inbam = new File(line.getOptionValue("i"));
		File outbam = new File(line.getOptionValue("o"));
		Integer orgQual = Integer.parseInt(line.getOptionValue("q"));
		Integer newQual = Integer.parseInt(line.getOptionValue("n"));
		adjustQualities(inbam, outbam, orgQual, newQual);
	}
	
	public static void adjustQualities(File infile, File outfile, int inputQual, int outputQual) throws IOException{
		SamReader reader = SamReaderFactory.makeDefault().open(infile);
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		SAMFileWriter writer = factory.makeBAMWriter(reader.getFileHeader(), true, outfile);
		int qualOffset = inputQual - outputQual;
		for (SAMRecord r : reader){
			byte[] orgQualities = r.getBaseQualities();
			byte[] newQualities = new byte[orgQualities.length];
			for (int i = 0; i < orgQualities.length; i++){
				newQualities[i] = (byte) (orgQualities[i] - qualOffset);
			}
			r.setBaseQualities(newQualities);
			writer.addAlignment(r);
		}
		reader.close();
		writer.close();
	}

}
