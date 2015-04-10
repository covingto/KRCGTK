package org.bcm.hgsc.cancer.pacbio;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
import org.bcm.hgsc.cancer.sv.Junction;
import org.bcm.hgsc.cancer.utils.Chromosome;
import org.bcm.hgsc.cancer.utils.Orientation;
import org.bcm.hgsc.utils.RangUtils;
import org.bcm.hgsc.utils.ReadFeature;
import org.bcm.hgsc.utils.Settings;

public class BreakMaker {

	public static int buffer = 100;
	private static int bedbuffer = 4000;
	/**
	 * @param args
	 * @throws IOException 
	 * @throws ParseException 
	 */
	public static void main(String[] args) throws IOException, ParseException {
		// TODO Auto-generated method stub
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("d", false, "Should debug be enabled");
		options.addOption("buffer", true, "How much buffer to use for including meta junctions [100]");
		options.addOption("bed", true, "Optional bed file of regions to consider, only breaks where one region falls in this bed file will be used.");
		options.addOption("bedbuffer", true, "If bed file is supplied how far away from a bed region can the break be to be considered [4000]");
		// options.addOption("germline", true, "If germline junctions are available, split these out into the germline file.");
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line = parser.parse(options, args);
		Settings.debug = line.hasOption("d");
		if (Settings.debug){ System.out.println(args); }
		String bedfilename = line.getOptionValue("bed", null);
		// String germlineblast = line.getOptionValue("germline", null);
		
		File bedfile = null;
		if (bedfilename != null){
			bedfile = new File(bedfilename);
			if (!bedfile.canRead()){
				System.out.println("Can't read " + bedfilename);
				System.exit(10);
			}
		}
		buffer = Integer.parseInt(line.getOptionValue("buffer", "100"));
		bedbuffer  = Integer.parseInt(line.getOptionValue("bedbuffer", "4000"));
		String[] reqargs = line.getArgs();
		
		if (reqargs.length < 2) {
			System.out.println("Required args not supplied.");
			formatter.printHelp("BreakMaker.jar", options);
			System.exit(10);
		}
		System.out.println("Starting analysis");
		System.out.println(Settings.dateFormat.format(new Date()));
		// File reffasta = new File(reqargs[0]);
		File splitreads = new File(reqargs[0]);
		String outbase = reqargs[1];
		
		// IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(reffasta);
		
		// build a library of reads from the blast outputs.  This library will then be filtered to identify junctions
		
		
		System.out.println("Building read library");
		// build a read library representing the reads and their mapping locations.
		System.out.println(Settings.dateFormat.format(new Date()));
		BlastStruct reads = new BlastStruct(splitreads);
		
		System.out.println("Blast result summary");
		reads.printSummary();
		
		System.out.println("Finding fully mapped reads");
		BlastGroup.filterFullyMapped(reads, new File(outbase + ".fullymapped"));
		reads.printSummary();
		
		System.out.println("Finding broken reads");
		// filter reads that are considered unmapped
		BlastGroup.filterUnmapped(reads, new File(outbase + ".unmapped"));
		
		System.out.println("Blast records post filter.");
		reads.printSummary();
		
		if (bedfile != null){
			System.out.println("Filtering using BED file.");
			BlastGroup.filterBed(reads, bedfile, bedbuffer);
		}
		
		// from the passing reads (those that might have junctions) we then make a junctions structure.
		// The junctions are then filtered and processed to make metajunctions (collections of multiple junctions).
		
		System.out.println("Building junctions");
		System.out.println(Settings.dateFormat.format(new Date()));
		List<Junction> junctionsa = buildJunctions(reads);
		
		
		
		System.out.println("Identifying recurrent breakpoints");
		List<Junction> recurrent = findRecurrentFusions(junctionsa, buffer);
		
		List<Junction> recurrents = recurrent;
		/*
		if (germlineblast != null){
			// This section handles the germline subtraction, we build a list
			File germlinefile = new File(germlineblast);
			if (!germlinefile.canRead()){
				System.out.println("Error reading germline file!");
				System.exit(10);
			} else {
				BlastStruct greads = new BlastStruct(germlinefile, fasta);
				List<Junction> gjunctions = buildJunctions(greads);
				recurrents = filterGermline(recurrent, gjunctions, new File(outbase + ".germline"));
			}
		} else {
			// fallback is to just assume that all is somatic or at least that the variable shoud retain.
			recurrents = recurrent;
		}
		*/
		
		System.out.println("Identified " + recurrents.size() + " recurrent fusions.");
		
		// print junctions to file
		System.out.println("Writing junctions.");
		System.out.println(Settings.dateFormat.format(new Date()));
		
		writeJunctions(junctionsa, new File(outbase + "_AllJunctions.txt"));
		writeJunctions(recurrents, new File(outbase + "_Junctions.txt"));
		
		
		// make metajunctions and score for support
		System.out.println("Building metajunctions.");
		System.out.println(Settings.dateFormat.format(new Date()));
		List<Metajunction> metajuctions = makeMetaJunctions(recurrents, buffer);
		
		// print junctions to file
		System.out.println("Writing metajunctions.");
		File metajunctionOuptutfile = new File(outbase + "_Metajunctions.txt");

		FileWriter fm = new FileWriter(metajunctionOuptutfile);
		BufferedWriter mjw = new BufferedWriter(fm);

		mjw.write("#support\tchra\tposa\tchrb\tposb\tjtype\trefgap\treadgap\torientation\treads\tastatus\tbstatus\n");
		for (Metajunction mj : metajuctions){
			mjw.write(mj.toString());
			mjw.newLine();
		}
		mjw.close();
		
		System.out.println("Building gap output.");
		// find gapped reads
		List<Junction> gaps = findGappedJunctions(recurrents);

		// print gapped junctions to file
		System.out.println("Writing gapped junctions.");
		File gapsOutputFile = new File(outbase + "_Gapped_Junctions.txt");

		FileWriter fg = new FileWriter(gapsOutputFile);
		BufferedWriter fgw = new BufferedWriter(fg);
		for (Junction gj : gaps){
			fgw.write(gj.juctionOuptut() + "\t" + gj.readGap() + "\t" + gj.referenceGap() + "\t" + gj.gapRatio() + "\t" + gj.queryPosLeft() + "\t" + gj.queryPosRight() + "\t" + gj.refPosLeft() + "\t" + gj.refPosRight());
			fgw.newLine();
		}
		fgw.close();
		
		// print seg dups
		System.out.println("Building tandem dups.");
		List<Junction> tanddups = findSegDups(recurrents);
		
		File tanddupOutputFile = new File(outbase + "_TandemDups_Junctions.txt");

		FileWriter fs = new FileWriter(tanddupOutputFile);
		BufferedWriter fsw = new BufferedWriter(fs);

		for (Junction gj : tanddups){
			fsw.write(gj.juctionOuptut() + "\t" + gj.readGap() + "\t" + gj.referenceGap() + "\t" + gj.gapRatio() + "\t" + gj.queryPosLeft() + "\t" + gj.queryPosRight() + "\t" + gj.refPosLeft() + "\t" + gj.refPosRight());
			fsw.newLine();
		}
		fsw.close();

		System.out.println("BreakMaker done.  Have a nice day!");
		System.out.println(Settings.dateFormat.format(new Date()));
	}
	
	public static void writeJunctions(List<Junction> junctions, File outfile) throws IOException{
		
		FileWriter fw = new FileWriter(outfile);
		BufferedWriter writer = new BufferedWriter(fw);
		
		for (Junction j : junctions){
			writer.write(j.juctionOuptut());
			writer.newLine();
		}
		writer.close();
	}
	
	/*
	 * Splits out germline junctions from recurrent junction list.
	 * @param Candidate list of recurrent junctions.
	 * @param Junctions to filter against.
	 * @param Output file for filtered junctions.
	 * @return
	 * @throws IOException 
	 *
	private static List<Junction> filterGermline(List<Junction> recurrent,
			List<Junction> gjunctions, File file) throws IOException {
			FileWriter fw = new FileWriter(file);
			BufferedWriter fsw = new BufferedWriter(fw);
			for ( Junction j : recurrent ){
				
			}
		return null;
	}
	*/
	
	/**
	 * Function returns reads that contain tandem duplications.
	 * 		ref		=====================================================
	 * 		map1		----------------------
	 * 		map2						-------------------------
	 * 		or
	 * 		map2						------
	 * 		map3						-------------------------
	 * Tandem duplications in our case have the following properties;
	 * 		1. map1 qend and map2 qstart are very close together (+/- 10bp)
	 * 		2. map2 sstart is less than map1 send
	 * 		3. map2 send should either be after map1 send or should be very near map1 send (-10bp)
	 * 		4. map1 and map2 must be in the same orientation
	 * @param recurrent
	 * @return
	 */

	private static List<Junction> findSegDups(List<Junction> recurrent){
		// the recurrent list has only junctions that are recurrent, but they are not ordered by the read groups
		// first restructure into read groups 
		//Map<String, SortedSet<Junction>> groupedJunctions = new HashMap<String, SortedSet<Junction>>();
		List<Junction> jlist = new ArrayList<Junction>();
		
		// collect all of the recurrent 
		for ( Junction j : recurrent ){
			if (!j.a.chr.equals(j.b.chr)){ continue; } // fails assertion 2
			if (!j.ol.equals(j.or)){ continue; } // fails assertion 4
			if (RangUtils.between(0.1, 10, Math.abs(j.gapRatio()))){ continue; } // fails assertion 1
			if (Math.abs(j.brl.qend - j.brr.qstart) > 10){ continue; }
			if (j.brl.send < j.brr.sstart){ continue; } // fails assertion 2
			if (j.brl.send - j.brr.send > 10){ continue; } // fails assertion 3
			jlist.add(j);
			//String brl = j.brl.q;
			// add a new entry for the query if it does not exist, sorting by start location of the query
			//if ( ! groupedJunctions.containsKey(brl) ){
			//	groupedJunctions.put(brl, new TreeSet<Junction>(Junction.qstartcompare));
			//}
			// exclude regions that are not sufficiently gapped, proper tandem dups will represent large differences in size between the query (which should be perfect) and the subject.
			//if (! RangUtils.between(0.1, 10, Math.abs(j.gapRatio()))){
			//	groupedJunctions.get(brl).add(j);
			//}
			
		}
		
		// now parse throught that structure and find candidate junctions
		/*for ( String k : groupedJunctions.keySet() ){
			SortedSet<Junction> jset = groupedJunctions.get(k);
			for (Junction j : jset){
				
			}
			
			Junction jf = jset.first();
			Junction jl = jset.last();
			BlastRow bf = jf.invert ? jf.brr : jf.brl;
			BlastRow bl = jl.invert ? jl.brr : jf.brl;
			if (bf.schr.equals(bl.schr)){
				// we got one, the chr are the same, this is a candidate insertion
			}
			
			jlist.addAll(jset);
			
		}*/
		return jlist;
	}

	private static List<Junction> findGappedJunctions(List<Junction> recurrent) {
		List<Junction> gappedJunctions = new ArrayList<Junction>();
		for (int i = 0; i < recurrent.size(); i++){
			Junction ju = recurrent.get(i);
			// potential gap if intrachrom and orientations match
			// other combinations are automatic SV candidates, not indel candidates, longer
			if ((ju.fusionType == Junction.FusionType.INTRACHROM) && (ju.ol == ju.or)){
				gappedJunctions.add(ju);
			}
		}
		return gappedJunctions;
	}

	private static List<ReadFeature> findRecurrentTerminals(List<ReadFeature> breaksBed, int buffer2){
		
		Set<ReadFeature> result		= new TreeSet<ReadFeature>();

		if (Settings.debug){
			System.out.println("Searching breaks");
		}

		for (int i = 0; i < breaksBed.size() - 1; i++ ){
			if ( Settings.debug & i % 1000 == 0){
				System.out.print("Processed " + i + " of " + breaksBed.size() + " total breaks\r");
			}
			int j = i + 1;
			ReadFeature a = breaksBed.get(i);
			while( j < breaksBed.size() ){
				ReadFeature b = breaksBed.get(j);
				
				if (! a.chr.equals(b.chr) ){ break; } // reached the end of the chromosome, no more matches possible
				if ( b.getStart() - a.getStart() > buffer2 ){ break; } // reached the end of the buffer, no more matches possible
				if (! a.read.equals(b.read) ){
					// we would be interested in adding these
					result.add(a);
					result.add(b);
					if (Settings.debug){
						System.out.print(a.toString() + ">" + b.toString() + "\r");
					}
					break;
				}
				j++;
			}
		}
		
		List<ReadFeature> recurrent = new ArrayList<ReadFeature>(result.size());
		recurrent.addAll(result);
		return recurrent;
	}
	
	private static List<ReadFeature> buildReadFeatures(List<Junction> junctions){
		// build a set of breaks
		TreeSet<ReadFeature> breaksSet 	= new TreeSet<ReadFeature>();

		if (Settings.debug){
			System.out.println("Building breaks Bed file");
		}

		for (Junction j : junctions){
			breaksSet.add(j.a);
			breaksSet.add(j.b);
		}
		List<ReadFeature> breaksBed = new ArrayList<ReadFeature>(breaksSet.size());
		breaksBed.addAll(breaksSet);
		return breaksBed;
	}
	
	/**
	 * Loop through all of the junctions and find those that are near to each other.  Only junctions supported by two events are returned.
	 * @param junctions
	 * @param buffer
	 * @return junctions List<Junction>
	 */
	public static List<Junction> findRecurrentFusions(List<Junction> junctions, int buffer){
		Set<Junction> resultSet = new TreeSet<Junction>();
		
		for (int i = 0; i < junctions.size() - 1; i++){
			Junction a = junctions.get(i);
			if (Settings.debug && i % 1000 == 0){
				System.out.print("Processed " + i + " junctions.\r");
				if (i % 100000 == 0){ System.out.println("Processed " + i + " junctions."); }
			}
			for (int j = i + 1; j < junctions.size(); j++){
				Junction b = junctions.get(j);
				if (a.zmw.equals(b.zmw)){ continue; } // keep itterating we aren't out of this read group yet
				if (a.a.near(b.a, buffer)){
					if (a.like(b, buffer)){
						resultSet.add(a);
						resultSet.add(b);
					}
				} else {
					break;
				}
			}
		}
		
		List<Junction> resultList = new ArrayList<Junction>(resultSet.size());
		resultList.addAll(resultSet);
		
		return resultList;
	}
	
	public static List<Junction> findRecurrentBreakpoints(List<Junction> junctions,
			int buffer2) {
		
		// generate a list of read features
		List<ReadFeature> breaksBed = buildReadFeatures(junctions);
		if (Settings.debug){
			System.out.println("Built " + breaksBed.size() + " breaks");
		}
		
		List<ReadFeature> recurrent = findRecurrentTerminals(breaksBed, buffer2);
		if (Settings.debug){
			System.out.println("Built " + breaksBed.size() + " breaks");
		}
		// 
		
		return filterJunctions(junctions, recurrent);
	}
	
//	private static boolean reverseContainsReadFeature(List<? extends BEDRegion> regions, BEDRegion r){
//		for (int i = regions.size() - 1; i >= 0; i--){
//			BEDRegion a = regions.get(i);
//			if (a.equals(r)){ return true; }
//			if (!a.getSequence().equals(r.getSequence())){ return false; }
//		}
//		return false;
//	}

	private static List<Junction> filterJunctions(List<Junction> junctions,
			List<ReadFeature> recurrent) {
		
		// build a set of hashes for the juntions to be added to, these can be grown efficiently
		Set<Junction> passingJunctions = new HashSet<Junction>();
		Set<Junction> junctionscopy = new HashSet<Junction>();
		junctionscopy.addAll(junctions); // now it's as a hashset
		
		for (int i = 0; i < recurrent.size(); i++){
			ReadFeature rf = recurrent.get(i);
			if ( Settings.debug && i % 1000 == 0 ){
				System.out.print("Processed " + i + " of " + recurrent.size() + " with " + passingJunctions.size() + " passing junctions");
			}
			for (Junction j : junctionscopy){
				if (j.a.equals(rf) || j.b.equals(rf)){
					passingJunctions.add(j);
					junctionscopy.remove(j);
					break;
				}
			}
		}
		
		List<Junction> passingList = new ArrayList<Junction>(passingJunctions.size());
		passingList.addAll(passingJunctions);
		return passingList;
	}

	private static List<Metajunction> makeMetaJunctions(List<Junction> junctions, int b) {
		List<Metajunction> mj = new ArrayList<Metajunction>();
		Set<Junction> seen = new TreeSet<Junction>();
		
		
		// need to put junctions into metajunctions
		for (int i = 0; i < junctions.size() - 1; i++){
			Junction seed = junctions.get(i);
			if (seen.contains(seed)){ continue; } // can't seed something that is already in a metajunction
			Metajunction candidate = new Metajunction(seed, b);
			for (int j = i + 1; j < junctions.size(); j++){
				Junction test = junctions.get(j);
				if ( seed.a.read.equals(test.a.read) ){ continue; } // can't make metajunction from same read
				if ( ! seed.a.chr.equals(test.a.chr) ){ break; } // no more metajunctions possible
				if ( candidate.fits(test) ){
					seen.add(test);
					candidate.add(test);
				}
			}
			if (candidate.size() > 1){
				mj.add(candidate);
			}
			
		}
		
		return mj;
	}

	private static List<Junction> buildJunctions(BlastStruct reads) {
		
		List<Junction> junctions = new ArrayList<Junction>();
		
		// loop across the reads
		System.out.println("Constructing junctions.");
		System.out.println(Settings.dateFormat.format(new Date()));
		for (String k : reads.getMappings().keySet()) {
			BlastGroup bg = reads.getMappings().get(k);
			if (bg.mappingLength() < 2){ continue; } // no need to process if only 1 read, findJunctions function should not be bothered by this but continue just in case
			junctions.addAll(findJunctions(bg));
		}
		
		// sort putative junctions
		System.out.println("Sorting junctions.");
		System.out.println(Settings.dateFormat.format(new Date()));
		Collections.sort(junctions, Junction.jcomparable);
		
		return junctions;
	}
	
	private static int maxclip(BlastRow rb, double clippct){
		// calculate the size of this mapping
		int dist = rb.qend - rb.qstart;
		
		// return the position within the read representing the start of the mapping + clippct of the size
		//             |------------|
		// ------------------------------------
		//                         *->
		return (int) Math.round(rb.qstart + dist * clippct);
	}

	private static Collection<Junction> findJunctions(BlastGroup bg) {
		List<Junction> js = new ArrayList<Junction>();
		
		// start the maxclip at maxinteger, we will bring this down soon
		int maxclippos = Integer.MAX_VALUE;
		
		// if there aren't more than two mappings return the empty string
		if (bg.mappingLength() < 2) { return js; }
		
		// start building the junctions
		for (int i = 0; i < bg.mappingLength()-1; i++){ 
			BlastRow ra = bg.getMapping(i);
			
			// new mappings can not start before this point
			int minstart = maxclip(ra, 0.9);
			for (int j = i+1; j < bg.mappingLength(); j++){
				BlastRow rb = bg.getMapping(j);
				boolean invert = false;
				// we can't start making junctions yet
				// check that the junctions are not 90% overlap and that the positions are at least 200 bp away from eachother (start to start).
				if (rb.qstart < minstart | Math.abs(rb.qstart - ra.qstart) < 200){ continue; }
				if (rb.qstart > maxclippos){ break; }
				int testmaxclippos = maxclip(rb, 0.9);
				if (testmaxclippos < maxclippos) { maxclippos = testmaxclippos; }
				// we don't need to process if either read is redundant
				if (ra.unique == BlastRow.REDUNDANT || rb.unique == BlastRow.REDUNDANT){ continue; }
				// we don't need to process if the reads are entirely ambiguous of eachother (two different mappings of the same region).
				if (BlastRow.unique(ra, rb) == BlastRow.AMBIGUOUS){ continue; }
				BlastRow ja = null;
				BlastRow jb = null;
				// now check if we need to reorient the rows
				// case, both are positiive; keep the current orientation
				if (ra.o == Orientation.POS & rb.o == Orientation.POS){
					ja = ra;
					jb = rb;
				// case that both are negative, 
				} else if (ra.o == Orientation.NEG & rb.o == Orientation.NEG){
					ja = rb.invert();
					jb = ra.invert();
					invert = true;
				// case where discordant (other cases caught above), now check genomic position, lowest on the left
				} else {
					if (ra.isGreaterThan(rb)){
						ja = rb.invert();
						jb = ra.invert();
						invert = true;
					} else {
						ja = ra;
						jb = rb;
					}
				} 
				if (ja == null || jb == null){
					System.out.println("Error; could not reconsile orientations of;");
					System.out.println(ra.toString());
					System.out.println(rb.toString());
					continue;
				}
				Chromosome achr = new Chromosome(ja.schr);
				Chromosome bchr = new Chromosome(jb.schr);
				if (achr.getChrnum() < 1 || bchr.getChrnum() < 1){ continue; }
				js.add(new Junction(new Chromosome(ja.schr), new Chromosome(jb.schr), ja.send, jb.sstart, ja.o, jb.o, ja, jb, invert));
			}
		}
		return js;
	}
	

}
