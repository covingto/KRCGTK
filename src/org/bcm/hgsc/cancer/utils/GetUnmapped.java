package org.bcm.hgsc.cancer.utils;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
import org.bcm.hgsc.cancer.bed.BEDRegion;
import org.bcm.hgsc.cancer.bed.BedTools;

public class GetUnmapped {

	/**
	 * @param args
	 * @throws ParseException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws ParseException, IOException {
		// TODO Auto-generated method stub
		// command line options
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("bed", true, "Optional bed file to search the sam/bam with.  If supplied, unmapped reads will be extracted using the indicated bedfile");
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line = parser.parse(options, args);
		String bedfilepath = line.getOptionValue("bed", null);
		String[] reqargs = line.getArgs();
		
		if (reqargs.length < 1){
			System.out.println("Required args not supplied.");
			formatter.printHelp("CollectBlocks.jar", options);
			System.exit(10);
		}
		String bamfilepath = reqargs[0];
		
		findUnmapped(bamfilepath, bedfilepath);
	}

	private static void findUnmapped(String bamfilepath,
			String bedfilepath) throws IOException {
		File bamfile = new File(bamfilepath);
		if (!bamfile.canRead()) {
			throw new IOException("Can't read bam file");
		}
		if (bedfilepath != null){
			File bedfile = new File(bedfilepath);
			if (!bedfile.canRead()) { 
				throw new IOException("Can't read bed file");
			}
			findUnmappedBed(bamfile, bedfile);
		} else {
			findUnmappedNoBed(bamfile);
		}
		
		
	}

	private static void findUnmappedNoBed(File bamfilepath) throws IOException {
		// just read through the entire bam file, find the unmapped reads and output both records
		final SamReader inputSam = SamReaderFactory.makeDefault().open(bamfilepath);
		final SamReader supportSam = SamReaderFactory.makeDefault().open(bamfilepath);
		final SAMFileWriter outsam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
				true, new File(bamfilepath + ".unmapped.bam"));
		SAMRecordIterator sri = inputSam.iterator();
		
		while (sri.hasNext()){
			final SAMRecord sr = sri.next();
			if (sr.getReadUnmappedFlag()){ continue; }
			if (!sr.getMateUnmappedFlag()){ continue; }
			final SAMRecord mate = supportSam.queryMate(sr);
			outsam.addAlignment(sr);
			outsam.addAlignment(mate);
		}
		
		// tidy up the files
		outsam.close();
		inputSam.close();
		supportSam.close();
		
	}

	private static void findUnmappedBed(File bamfilepath, File bedfile) throws IOException {
		// TODO Auto-generated method stub
		List<BEDRegion> regions = BedTools.processBedRegions(null, bedfile.getAbsolutePath());
		final SamReader inputSam = SamReaderFactory.makeDefault().open(bamfilepath);
		final SamReader supportSam = SamReaderFactory.makeDefault().open(bamfilepath);
		final SAMFileWriter outsam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
				true, new File(bamfilepath + ".unmapped.bam"));
		for (BEDRegion br : regions){
			System.out.println("Processing region " + br.getSequence() + ":" + br.getStart() + "-" + br.getEnd());
			SAMRecordIterator sri = inputSam.query(br.getSequence(), br.getStart(), br.getEnd(), false);
			while (sri.hasNext()){
				final SAMRecord sr = sri.next();
				final SAMRecord mate = supportSam.queryMate(sr);
				if (mate.getMateUnmappedFlag()){
					outsam.addAlignment(sr);
					outsam.addAlignment(mate);
					
				}
			}
			sri.close();
		}
		
		// tidy up the files
		outsam.close();
		inputSam.close();
		supportSam.close();
				
	}

}
