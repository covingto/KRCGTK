package org.bcm.hgsc.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
import org.bcm.hgsc.cancer.bed.BEDRegion;
import org.bcm.hgsc.cancer.bed.BedTools;

public class FindPairs {
	public static int chunkbuffer = 1000;
	/**
	 * @param args
	 * @throws ParseException 
	 */
	public static void main(String[] args) {
		// parse command line arguments
		Options options = new Options();
		Parser parser = new BasicParser();
		// options.addOption("bed", true, "Optional bed file to search the sam/bam with.  If supplied, blocks will only be identified in the indicated bedfile");
		options.addOption("v", false, "Verbose output used");
		// options.addOption("chunkbuffer", true, "Local buffer to search for mates to the initial query, this speeds the lookup of the mates in the bam. [1000]");
		options.addOption("h", false, "Help");
		// options.addOption("h", false, "Should the bam header be appended [false]");
		HelpFormatter formatter = new HelpFormatter();
		String bedfilepath = null;
		String bamfile = null;
		try {
			CommandLine line = parser.parse(options, args);
			if (line.hasOption("h")){
				System.err.println(
						"FindPairs; \n Author; Kyle R. Covington (covingto@bcm.edu) \n" +
						" Purpose; Extracts pairs of reads where any read intersects with regions in the defined bed file.\n" +
						"   Reads are writen in SAM format (with no header) to files labeled reads<chr>_<posa>_<posb>.fp.sam\n" +
						"   where chr posa and posb are from the bed file regions."
						);
				System.exit(1);
			}
			String[] reqargs = line.getArgs();
			if (reqargs.length < 2){
				System.err.println("Required args not supplied.");
				formatter.printHelp("FindPairs.jar", options);
				System.exit(10);
			}
			bedfilepath = reqargs[0];
			bamfile = reqargs[1];
			FindPairs.chunkbuffer = Integer.parseInt(line.getOptionValue("chunkbuffer", "1000"));
			Settings.debug = line.hasOption("v");
		} catch (ParseException pe) {
			System.err.println("Error parsing arguments.");
			formatter.printHelp("FindPairs.jar", options);
			System.exit(10);
		}
		
		
		
		try {
			// start the processing of the bam file
			List<BEDRegion> regions = BedTools.processBedRegions(null, bedfilepath);
			if (Settings.debug){ System.out.println("Processing Bed regions"); }
			for (BEDRegion br : regions){
				List<PairedReads> pairs = _getLocalPairedReads(br, bamfile);
				BufferedWriter writer = new BufferedWriter(new FileWriter("reads" + br.getSequence() + "_" + br.getStart() + "_" + br.getEnd() + ".fp.sam")); 
				for (PairedReads p : pairs){
					p.write(writer);
				}
				writer.close();
				if (Settings.debug){ System.out.println("Wrote file; " + "reads" + br.getSequence() + "_" + br.getStart() + "_" + br.getEnd() + ".fp.sam"); }
			}
		} catch (IOException e) {
			System.err.println("IOException, please ensure that bedfile and bamfile are readable and outfile directory is writeable.");
			System.exit(10);
		}
		
	}
	
	public static List<PairedReads> getPairedReads(String bedfilepath, String bamfile) throws IOException{
		return getPairedReads(bedfilepath, bamfile, 1);
	}

	public static List<PairedReads> getPairedReads(String bedfilepath,
			String bamfile, int threads) throws IOException {
		List<PairedReads> pairs = new ArrayList<PairedReads>();
		List<BEDRegion> regions = BedTools.processBedRegions(null, bedfilepath);
		for (BEDRegion r : regions){
			pairs.addAll(_getLocalPairedReads(r, bamfile));
		}
		return null;
	}

	private static List<PairedReads> _getLocalPairedReads(
			BEDRegion r, String bamfile) throws IOException {
		List<PairedReads> pairs = new ArrayList<PairedReads>();
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		SamReader inputSam = SamReaderFactory.makeDefault().open(new File(bamfile));
		SAMRecordIterator sri = inputSam.query(r.getSequence(), r.getStart(), r.getEnd(), false);
		while (sri.hasNext()){
			final SAMRecord sr = sri.next();
			reads.add(sr);
		}
		sri.close();
		
		// pair up the mates
		boolean[] flagged = new boolean[reads.size()];
		Arrays.fill(flagged, false);
		for (int i = 0; i < reads.size() - 1; i++){
			if (flagged[i]) continue; // flagged means that we have already put the read into a pair and don't need to do that again
			final String readName = reads.get(i).getReadName();
			boolean added = false;
			for (int j = i + 1; j < reads.size(); j++){
				if (flagged[j]) continue;
				if (reads.get(j).getReadName().equals(readName)){
					pairs.add(new PairedReads(reads.get(i), reads.get(j)));
					flagged[i] = true;
					flagged[j] = true;
					added = true;
					break; // found the mate, stop looking for more
				}
			}
			// if the pair wasn't created in the section we looked up then check the entire bam for the mate
			if ( !added ){
				pairs.add(new PairedReads(reads.get(i), inputSam.queryMate(reads.get(i))));
				flagged[i] = true;
			}
		}
		inputSam.close();
		return pairs;
	}
	
	

}
