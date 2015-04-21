package org.bcm.hgsc.utils;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import net.minidev.json.JSONObject;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.bcm.hgsc.utils.BAMUtils.ConformedRead;

public class BAMQC {
	private final Map<String, Integer> cigarMap1 = new HashMap<String, Integer>();
	private final Map<String, Integer> cigarMap2 = new HashMap<String, Integer>();
	private final Map<String, Integer> varMap1 = new HashMap<String, Integer>();
	private final Map<String, Integer> varMap2 = new HashMap<String, Integer>();
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("bam", true, "BAM file");
		options.addOption("ref", true, "Reference fasta files, requires indexing");
		options.addOption("out", true, "Output file (json)");
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line = parser.parse(options, args);
		
		File bam = new File(line.getOptionValue("bam"));
		File output = new File(line.getOptionValue("out"));
		File ref = new File(line.getOptionValue("ref"));
		BAMQC bqc = new BAMQC();
		bqc.run(bam, ref, output);
	}
	
	public void run(File bam, File ref, File output) throws Exception{
		final int buffer = 1000000;
		// activate file readers
		final IndexedFastaSequenceFile fastaref = new IndexedFastaSequenceFile(ref);
		final SAMSequenceDictionary seq_dict = fastaref.getSequenceDictionary();
		final BAMInterface bam_interface = new BAMInterface(bam, bam.getName(), "test");
		
		// start the pool
		final ExecutorService pool = Executors.newFixedThreadPool(12);
		for (final SAMSequenceRecord r : seq_dict.getSequences()){
			int end_pos = 0;
			while ((end_pos + buffer) < r.getSequenceLength()){
				Runnable worker = new QCWorker(bam_interface, r.getSequenceName(), 
						end_pos, end_pos + buffer, this, fastaref);
				pool.execute(worker);
				end_pos = end_pos + buffer;
			}
			Runnable worker = new QCWorker(bam_interface, r.getSequenceName(), 
					end_pos, r.getSequenceLength(), this, fastaref);
			pool.execute(worker);
		}
		pool.shutdown();
		while (!pool.isShutdown()){
			Thread.sleep(10000);
		}
		
		this.printMaps(output);
	}
	
	public void printMaps(File output){
		JSONObject value = new JSONObject();
		value.put("r1_cigar", cigarMap1);
		value.put("r2_cigar", cigarMap2);
		value.put("r1_vars", varMap1);
		value.put("r2_vars", varMap2);
		try {
			PrintWriter out = new PrintWriter(output.getAbsolutePath());
			out.write(value.toJSONString());
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public class QCWorker implements Runnable {
		final BAMInterface bam_interface;
		final String sequenceName;
		final int start;
		final int end;
		final IndexedFastaSequenceFile fastaref;
		public boolean success = false;
		private final BAMQC bqc;
		private final Map<String, Integer> cigarMap1 = new HashMap<String, Integer>();
		private final Map<String, Integer> cigarMap2 = new HashMap<String, Integer>();
		private final Map<String, Integer> varMap1 = new HashMap<String, Integer>();
		private final Map<String, Integer> varMap2 = new HashMap<String, Integer>();
		public QCWorker(BAMInterface bam_interface, String sequenceName,
				int start, int end, BAMQC bqc, IndexedFastaSequenceFile fastaref) {
			this.bam_interface = bam_interface;
			this.sequenceName = sequenceName;
			this.start = start;
			this.end = end;
			this.bqc = bqc;
			this.fastaref = fastaref;
		}
		
		@Override
		public void run() {
			SamReader sam_reader = bam_interface.getSamfilereader();
			SAMRecordIterator sri = sam_reader.query(this.sequenceName, this.start, this.end, false);
			try{
				while (sri.hasNext()){
					final SAMRecord sr = sri.next();
					if (sr.getAlignmentStart() < this.start || sr.getAlignmentEnd() >= this.end){ continue; } // not really in this set
					if (!passingRead(sr)){ continue; }
					// handle the CIGAR
					final String cigar = sr.getCigarString();
					addCigar(sr.getFirstOfPairFlag() ? cigarMap1 : cigarMap2, cigar);
					// find the variants
					ConformedRead cr = BAMUtils.conformToReference(sr, fastaref);
					addVars(sr.getFirstOfPairFlag() ? varMap1 : varMap2, cr);
				}
				updateBQCHashes();
				this.success = true;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} finally {
				sri.close();
				try {
					sam_reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
		}
		
		private void updateBQCHashes() {
			updateHash(bqc.cigarMap1, this.cigarMap1);
			updateHash(bqc.cigarMap2, this.cigarMap2);
			updateHash(bqc.varMap1, this.varMap1);
			updateHash(bqc.varMap2, this.varMap2);
		}

		private void updateHash(Map<String, Integer> pmap,
				Map<String, Integer> cmap) {
			synchronized(pmap){
				for (Entry<String, Integer> entry : cmap.entrySet()){
					if (!pmap.containsKey(entry.getKey())){
						pmap.put(entry.getKey(), entry.getValue());
					} else {
						pmap.put(entry.getKey(), pmap.get(entry.getKey()) + entry.getValue());
					}
				}
			}
			
		}

		private void addVars(Map<String, Integer> map, ConformedRead cr) {
			// scan the conformed read for variants
			int lmp = 0;
			boolean extending = false;
			for (int i = 0; i < cr.read.length; i++){
				final byte read = cr.read[i];
				final byte ref = cr.ref[i];
				final CigarOperator co = cr.cigar[i];
				if (extending){
					// have not yet opened a gap
					if (co.equals(CigarOperator.S)){ break; }
					else if (co.equals(CigarOperator.M) && read == ref){
						// break extension
						addVar(map, Arrays.copyOfRange(cr.read, lmp, i + 1), Arrays.copyOfRange(cr.ref, lmp, i + 1));
						lmp = i;
					}
				} else {
					if (co.equals(CigarOperator.S)){ continue; }
					else if (!(co.equals(CigarOperator.M) && read == ref)){
						extending = true;
					} else {
						lmp = i;
					}
				}
			}
		}

		private void addVar(Map<String, Integer> map, byte[] read, byte[] ref) {
			final String allele = new String(read) + ">" + new String(ref);
			synchronized(map){
				if (!map.containsKey(allele)){
					map.put(allele, 1);
				} else {
					map.put(allele, map.get(allele) + 1);
				}
			}
			
		}

		private boolean passingRead(SAMRecord sr) {
			if (sr.getDuplicateReadFlag()){ return false; }
			else if (sr.getReadUnmappedFlag()){ return false; }
			else if (sr.getNotPrimaryAlignmentFlag()){ return false; }
			else if (sr.getReadFailsVendorQualityCheckFlag()){ return false; }
			return true;
		}

		private void addCigar(Map<String, Integer> cigarMap, String cigar){
			synchronized(cigarMap){
				if (!cigarMap.containsKey(cigar)){
					cigarMap.put(cigar, 1);
				} else {
					cigarMap.put(cigar, cigarMap.get(cigar) + 1);
				}
			}
		}
		
	}

}
