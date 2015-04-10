package org.bcm.hgsc.cancer.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.cancer.utils.LimitingStreamingWalkAcceptor;
import org.bcm.hgsc.cancer.utils.StreamingWalkAcceptor;
import org.bcm.hgsc.cancer.utils.WalkAcceptor;
import org.bcm.hgsc.utils.BAMInterface;
import org.bcm.hgsc.utils.Settings;

public class RefWalker {
	int lastpos = 0;
	static final char[] bases = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'n', 'N'};
	static final char[] numbers = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	static final char[] mdredirects = {'^'};
	public static int min_reads = 50;
	boolean test = false;
	
	SAMRecord lastRecord = null;
	/**
	 * This utility walks the reference and generates a set of channels representing the level of variation
	 * in a sliding window along the genome.  Window sizes and spacings can be set by the user.
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("t", false, "Should we run in test mode?  Prints results after processing from 10M to 11M on the first chr.");
		options.addOption("d", false, "Should debug be enabled");
		options.addOption("w", true, "Window size [1000]");
		options.addOption("b", true, "Breaks per window. [10]");
		//options.addOption("s", false, "Stream data from the bedGraphs, this option allows the process to proceed while consuming less memory since values are written as they are being read.  Averages are reported at the end of processing.");
		options.addOption("max", true, "Max values to report in bedGraph");
		options.addOption("min", true, "Min values to report in bedGraph");
		options.addOption("min_reads", true, "Min number of reads to calculate a precentage");
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line = parser.parse(options, args);
		
		Settings.debug = line.hasOption("d");
		int windowSize = Integer.parseInt(line.getOptionValue("w", "1000"));
		int breaks = Integer.parseInt(line.getOptionValue("b", "10"));
		double minval = Double.parseDouble(line.getOptionValue("min", "0"));
		double maxval = Double.parseDouble(line.getOptionValue("max", "-1"));
		min_reads = Integer.parseInt(line.getOptionValue("min_reads", "50"));
		String[] reqargs = line.getArgs();
		
		if (reqargs.length < 2){
			System.out.println("Required args not supplied.");
			formatter.printHelp("CollectBlocks.jar", options);
			System.exit(10);
		}
		String bamfile = reqargs[0];
		String outfileext = reqargs[1];
		
		RefWalker rw = new RefWalker();
		rw.test = line.hasOption("t");
		//if (stream){
		StreamingWalkAcceptor signalMM = new LimitingStreamingWalkAcceptor(new File(outfileext + "_MM.bedg"), minval, maxval);
		StreamingWalkAcceptor signalIM = new LimitingStreamingWalkAcceptor(new File(outfileext + "_IM.bedg"), minval, maxval);
		StreamingWalkAcceptor signalT = new LimitingStreamingWalkAcceptor(new File(outfileext + "_T.bedg"), minval, maxval);
		rw.walkBAM(bamfile, outfileext, windowSize, breaks, signalMM, signalIM, signalT);
		signalIM.close();
		signalMM.close();
		signalT.close();
		//} else {
		//	CollectingWalkAcceptor signalMM = new CollectingWalkAcceptor(new File(outfileext + "_MM.bedg"));
	//		CollectingWalkAcceptor signalIM = new CollectingWalkAcceptor();
		//	CollectingWalkAcceptor signalT = new CollectingWalkAcceptor();
		//	rw.walkBAM(bamfile, outfileext, windowSize, breaks, signalMM, signalIM, signalT);
		//}
	}
	
	public void walkBAM(String bamfile, String outfileext, int windowSize, int breaks, WalkAcceptor signalMM, WalkAcceptor signalIM, WalkAcceptor signalT) throws Exception{
		SamReader sfr = new BAMInterface(new File(bamfile), null, null).getSamfilereader();
		int blockSize = windowSize / breaks;


		// get a listing of the sequence names
		List<SAMSequenceRecord> seqrec = sfr.getFileHeader().getSequenceDictionary().getSequences();
		for (SAMSequenceRecord sr : seqrec){
			this.lastpos = 0;
			this.lastRecord = null;
			// steamrollers
			int[] mmBlocks = new int[breaks];
			int[] imBlocks = new int[breaks];
			long[] positions = new long[breaks];
			int[] counts = new int[breaks];
			
			int breakpos = this.test ? 10000000 : 0;
			int posindex = breaks / 2;
			int updateindex = 0;
			int positionNumber = 0;
			final String seqname = sr.getSequenceName();
			if (Settings.debug){ System.err.println("Parsing sequence " + seqname); }
			int start = this.test ? 10000000 : 0;
			int end   = this.test ? 11000000 : 0;
			SAMRecordIterator sri = sfr.query(seqname, start, end, false);
			if (Settings.debug){ System.err.println("Itterator generated"); }
			// populate the steam-rollers
			for (int i = 0; i < breaks; i++){
				breakpos += blockSize;
				List<SAMRecord> records = getRecords(sri, breakpos);
				mmBlocks[i] = calcMM(records);
				imBlocks[i] = calcIM(records);
				counts[i] = records.size();
				positions[i] = breakpos;
			}
			if (Settings.debug){ System.err.println("Initial block filling done"); }
			do {
				positionNumber++;
				if (Settings.debug && positionNumber % 100000 == 0){ System.err.println("Processed " + positionNumber + " blocks"); }
				final long thispos = positions[posindex];
				

				// calculate the average mm
				int mm = 0;
				int im = 0;
				int count = 0;
				for (int j = 0; j < breaks; j++){
					mm += mmBlocks[j];
					im += imBlocks[j];
					count += counts[j];
				}
				
				if (count > min_reads){
					double mmval = (double) mm / ((double) count * 100);
					double imval = (double) im / (double) count;
					double totval = ((double) mm + (im * 100)) / ((double) count * 100);
				
					if (this.test){
						System.out.println(StringUtils.join(mmBlocks, ","));
						System.out.println(StringUtils.join(imBlocks, ","));
						System.out.println(StringUtils.join(counts, ","));
						System.out.println("MM: " + mmval + " IM: " + imval + " TOT: " + totval + " Count: " + count);
					}
					// update the data for the current position
					signalMM.add(seqname, thispos, thispos + blockSize, mmval);
					signalIM.add(seqname, thispos, thispos + blockSize, imval);
					signalT.add(seqname, thispos, thispos + blockSize, totval);
					
				}
				
				List<SAMRecord> records = getRecords(sri, breakpos);
				mmBlocks[updateindex] = calcMM(records);
				imBlocks[updateindex] = calcIM(records);
				counts[updateindex] = records.size();
				positions[updateindex] = breakpos;
				posindex++;
				if (posindex >= breaks){ posindex -= breaks; }
				updateindex++;
				if (updateindex >= breaks){ updateindex -= breaks; }
				breakpos += blockSize;
			} while(sri.hasNext());
			sri.close();
			if (this.test){ break; } 
		}
		sfr.close();
	}
	
	private int calcIM(List<SAMRecord> records){
		int im = 0;
		for (SAMRecord sr : records){
			if (sr.getMateUnmappedFlag()){
				im++;
			} else if (!sr.getProperPairFlag()){
				im++;
			}
		}
		return im;
	}
	
	/**
	 * Calculate the number of non-reference bases in the given reads.
	 * @param records
	 * @return int, number of Match differecne, these include mis-matches, indels, and soft clipped bases
	 */
	private int calcMM(List<SAMRecord> records){
		int mm = 0;
		for (SAMRecord sr : records){
			final int mis = countMisMatches(sr);
			if (mis < 0) continue; // a -1 return value indicates that no MD tag was present, we continue the records
			final int ids = countIDS(sr);
			mm += mis + ids;
		}
		return mm;
		
	}
	
	private int countIDS(SAMRecord sr){
		return countIDS(sr.getCigar());
	}
	
	private int countIDS(Cigar cigar){
		int count = 0;
		for (final CigarElement cigEl : cigar.getCigarElements()){
			final int cigElLen = cigEl.getLength();
			final CigarOperator cigElOp = cigEl.getOperator();
			if (cigElOp == CigarOperator.DELETION){
				count += cigElLen;
			} else if (cigElOp == CigarOperator.INSERTION){
				count += cigElLen;
			} else if (cigElOp == CigarOperator.SOFT_CLIP){
				count += cigElLen;
			}
		}
		return count;
	}
	
	private int countMisMatches(SAMRecord sr){
		return countMisMatches(sr.getStringAttribute(SAMTag.MD.name()));
	}
	
	private int countMisMatches(String md){
		if (md == null) return -1;
		char[] chars = md.toCharArray();
		int count = 0;
		boolean inredirect = false;
		for (char c : chars){
			if (inredirect){
				if (charin(c, RefWalker.numbers)){
					inredirect = false;
					continue;
				}
			} else {
				if (charin(c, RefWalker.numbers)) continue;
				if (charin(c, RefWalker.mdredirects)){
					inredirect = true;
					continue;
				}
				if (charin(c, RefWalker.bases)){
					count++;
					continue;
				}
			}
		}
		return count;
	}
	
	private boolean charin(char c, char[] test){
		for (char t : test){
			if (c == t) return true;
		}
		return false;
	}
	
	private List<SAMRecord> getRecords(SAMRecordIterator sri, int end){
		List<SAMRecord> records = new ArrayList<SAMRecord>();
		if (lastRecord != null && lastpos < end){
			records.add(lastRecord);
			lastRecord = null;
		}
		if (lastpos >= end){
			// this means that even the last itteration would not get us to this block,
			// the records are empty
			return records;
		}
		while(sri.hasNext() && lastpos < end){
			final SAMRecord sr = sri.next();
			if (sr.getDuplicateReadFlag()) continue;
			if (sr.getReadUnmappedFlag()) continue;
			if (sr.getMappingQuality() < 20) continue;
			lastpos = sr.getAlignmentEnd();
			records.add(sr);
		}
		if (Settings.debug){ System.out.println("Fetched " + records.size() + " records"); }
		return records;
	}

}
