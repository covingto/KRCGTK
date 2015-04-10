package org.bcm.hgsc.cancer.virus;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingDeque;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
import org.bcm.hgsc.utils.BAMInterface;
import org.bcm.hgsc.utils.Settings;
import org.json.simple.JSONValue;

import com.skjegstad.utils.BloomFilter;

public class ViralBloomFilter {
	public static BloomWorker[] bloomWorkers = new BloomWorker[4];
	public static LinkedBlockingDeque<SAMRecord> samRecords = new LinkedBlockingDeque<SAMRecord>(100);
	public static Set<String> writtenUnmappedMates = new HashSet<String>();
	public static Set<SAMRecord> matesToGet = new HashSet<SAMRecord>();
	public static int hitCount = 0;

	/**
	 * Runs in two different modes, train and filter
	 * Train (-t).
	 * 	Train takes a configuration file and begins training the filter.  Input configuration file should be a json formatted file containing a "hash" with values.
	 * 	<ul>
	 * 		<li>bases - number of bases in the fasta, this is used to calculate the size of the filter</li>
	 * 		<li>files - location of the files used for filtering</li>
	 * 		<li>error (optional) - accepted error rate for the filter</li>
	 * 		<li>output - location of the output of the training.</li>
	 * 	</ul>
	 * Filter (default).
	 * 	Filter takes a BAM file and the result of training (as specified in output above) and triages unmapped reads based on presence in the filter.
	 * By default a SAM file is written containing the unmapped but likely hit read and it's mate, though this behavior can be altered with options.
	 * Options.
	 * 	-v	Collect unmapped reads that do not hit the filter (and their mates).
	 *  -m	Do not collect the mates of the reads.
	 * @param args
	 * @throws Exception 
	 * @throws ParseException 
	 * @throws ClassNotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws Exception{
		// Accepts two arguments, filter and train each with different inputs 
		// train is the first step.  This will generate a bloom filter based on input configurations and train the filter based on that
		// filter is the second step.  
		try{
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("t", true, "Train the filter, overrides all other arguments and only trains");
		options.addOption("v", false, "Collect unmapped reads that do not hit the filter");
		options.addOption("d", false, "Debug mode");
		options.addOption("m", false, "Do not output the mates");
		options.addOption("test", true, "Tests the given string to see if it has a hit with the database.");
		options.addOption("h", false, "Prints help and exits");
		// options.addOption("germline", true, "If germline junctions are available, split these out into the germline file.");
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line = parser.parse(options, args);
		Settings.debug = line.hasOption("d");

		if (line.hasOption("t")){
			_train(line.getOptionValue("t"));
		} else if (line.hasOption("h")) {
			//_useage(formatter);
			formatter.printHelp("ViralBloomFilter.jar", options);
		} else if (line.hasOption("test")){
			String[] reqargs = line.getArgs();
			_test(reqargs[0], line.getOptionValue("test"));
		} else {
			String[] reqargs = line.getArgs();
			if (reqargs.length < 3){
			//	_useage(formatter);
				formatter.printHelp("ViralBloomFilter.jar", options);
			} else {
				_filter(reqargs[0], reqargs[1], reqargs[2]);
			}
		}
		} catch  (IOException e){
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ParseException e) {
			System.out.println("Parser not properly formatted");
			e.printStackTrace();
		}
	}
	
	private static void _test(String filterFile, String testString) throws IOException, ClassNotFoundException {
		BloomTester tester = new BloomTester(filterFile);
		if (tester.test(testString)){
			System.out.println(testString + "\tFound");
		} else {
			System.out.println(testString + "\tNotFound");
		}
		
	}

	public static void _useage(HelpFormatter formatter){
		// TODO: add
	}

	@SuppressWarnings("unchecked")
	public static void _train(String config) throws IOException{
		FileReader fr;
		fr = new FileReader(config);

		Map<String, Object> conf = (Map<String, Object>) JSONValue.parse(fr);
		fr.close();
		
		List<String> targetGenomes = (List<String>) conf.get("genomes");
		String outFile = (String) conf.get("outfile");
		Float fpr = Float.parseFloat((String) conf.get("false_positive_rate"));
		Integer bins = Integer.parseInt((String) conf.get("bins"));
		
		BloomFilter<String> filter = new BloomFilter<String>(fpr, bins);
		
		ExecutorService pool = Executors.newFixedThreadPool(8);
		
		int nContigs = 0;
		for (String gf : targetGenomes){
			nContigs++;
			final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(gf));
			ReferenceSequence refseq;
			while ((refseq = fasta.nextSequence()) != null){
				pool.execute(new BloomTrainer(refseq, filter));
			}
			fasta.close();
		}
		
		System.out.println(nContigs + " submitted to the job pool");
		
		pool.shutdown();
		while (!pool.isTerminated()){
			try {
				Thread.sleep(60);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		System.out.println("Writing bloom filter to file " + outFile);
		FileOutputStream fileOut = new FileOutputStream(outFile);
		ObjectOutputStream out = new ObjectOutputStream(fileOut);
		out.writeObject(filter);
		out.close();
		System.out.println("Training finished, have a nice day!");
		
	}

	public static void _filter(String filterFile, String bamFile, String outputFile) throws Exception{
		//BloomFilter<String> filter = _deserialize(filterFile);
		System.out.println("Filtering " + bamFile + " through the bloom filter " + filterFile + "\nOutput will be in " + outputFile);
		SamReader inputSam = new BAMInterface(new File(bamFile), null, null).getSamfilereader();
		SynchronizedSamWriter writer = new SynchronizedSamWriter(outputFile, inputSam);
		// initialize the BloomWorkers
		System.out.println("Building Bloom workers");
		int workers = 4;
		int lines = 0;
		ExecutorService pool = Executors.newFixedThreadPool(workers);
		for (int w = 0; w < workers; w++){
			pool.execute(new BloomWorker(filterFile, bamFile, writer));
		}
		// start working
		System.out.println("Iterating over BAM");
		SAMRecordIterator sri = inputSam.queryUnmapped();
		while (sri.hasNext()){
			final SAMRecord sr = sri.next();
			if (matesToGet.contains(sr)){
				writer.writeHit(sr);
				matesToGet.remove(sr);
				continue;
			} 
			lines++;
			if (lines % 1000 == 0){
				System.out.println("Processed " + lines + " unmapped reads and got " + hitCount + " hits");
			}
			samRecords.put(sr);
		}
		System.out.println("Processed " + lines + " unmapped reads and got " + hitCount + " hits. Shutting down the pool.");
		pool.shutdown();
		while (!pool.isTerminated()){
			try {
				Thread.sleep(60);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		System.out.println("Workers done, writing unprocessed mates");
		// write the as yet unwritten mates
		for (SAMRecord sr : matesToGet){
			writer.writeHit(sr);
		}
		writer.close();
		System.out.println("Done, have a nice day");
		// done
	}

	@SuppressWarnings("unchecked")
	public static BloomFilter<String> _deserialize(String filterFile) throws IOException, ClassNotFoundException{
		InputStream file = new FileInputStream(filterFile);
		ObjectInput input = new ObjectInputStream(file);
		BloomFilter<String> filter;
		try{
			filter = (BloomFilter<String>) input.readObject();
		} finally {
			input.close();
		}
		return filter;
	}
}
