package org.bcm.hgsc.cancer.sv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.bcm.hgsc.cancer.utils.Orientation;

public class HoneyTailsJunctionReader implements JunctionSourceReader {

	private final BufferedReader reader;
	private JRecord next;
	private final Caller caller;
	
	public HoneyTailsJunctionReader(String source, String params, int buffer, Caller caller) throws FileNotFoundException{
		this.caller = caller;
		this.reader = new BufferedReader(new FileReader(new File(source)));
		this.next = parseLine();
	}
	
	public HoneyTailsJunctionReader(String source, Caller caller) throws FileNotFoundException{
		this(source, "default", 1000, caller);
	}
	
	/**
	 * pie syntax.
	 * ->p=i->=->i=e->
	 * 
	 * This syntax was (developed?) used by A. English for Honey, and can be a bit maddening at times.
	 * Any read (in read space) can be written as ->p=i->=->i=e-> where i indicates the junction start or end position of the 
	 * initial mapped section of the read.  p indicates the prolog or a bit that was initially not mapped but was then
	 * mapped by Honey-tails, and the e indicates the epilogue.  Upon rearrangement, several things can happen;
	 * 
	 * Positive positive (or negative negative) junctions;
	 *    =========>>>>============>
	 *    Forward reads;
	 *            ->p=i->
	 *            ->i=e->
	 *    Reverse reads;
	 *            <-e=i<-
	 *            <-i=p<-
	 *            
	 * Head to head junctions;
	 *    =========>><<=============
	 *    Forward reads;
	 *            ->i%<-e
	 *            <-e%->i (don't ask me why)
	 *            
	 *    Reverse reads;
	 *            <-i%->p
	 *            ->p%<-i
	 * 
	 * Tail to tail junctions;
	 *    <========<<>>=============>
	 *    Forward reads;
	 *            i->%p<-
	 *            p<-%i->
	 *    Reverse reads;
	 *            i<-%e->
	 *            e->%i<-
	 *            
	 * @param s
	 * @return
	 */
	private Orientation[] parseOrientations(String s){
		boolean strandSplit = s.contains("%");
		boolean strandNoSplit = s.contains("=");
		if (strandSplit && strandNoSplit){
			// we both split and didn't split
			return new Orientation[] {Orientation.AMBIGUOUS, Orientation.AMBIGUOUS};
		} else if (strandNoSplit){
			// these are easier cases since a negative negative orientaiton implies a positive positive one
			return new Orientation[] {Orientation.POS, Orientation.POS};
		} else {
			// darn this is an orientation split, we need to find out if it's head to head or tail to tail or both
			boolean headToHead = s.contains("->i%<-e") || s.contains("<-e%->i") || s.contains("<-i%->p") || s.contains("->p%<-i");
			boolean tailToTail = s.contains("i->%p<-") || s.contains("p<-%i->") || s.contains("i<-%e->") || s.contains("e->%i<-");
			if (headToHead && tailToTail){
				// both head to head and tail to tail junction
				return new Orientation[] {Orientation.RECIP, Orientation.RECIP};
			} else if (headToHead){
				return new Orientation[] {Orientation.POS, Orientation.NEG};
			} else if (tailToTail){
				return new Orientation[] {Orientation.NEG, Orientation.POS};
			} else {
				System.err.println("Could not parse orientation information; " + s);
				return new Orientation[] {Orientation.AMBIGUOUS, Orientation.AMBIGUOUS};
			}
		}
	}
	
	
	private JRecord parseLine(){
		String line;
		try {
			line = reader.readLine();
			if (line == null){
				return null;
			} else if (line.startsWith("#")){
				return parseLine(); // read the next line, we skip lines with #
			} else {
				String[] lsplit = line.split("\t");
				if (lsplit.length < 12){
					System.err.println("Line did not have 12 fields");
					System.err.print(line);
					return parseLine();
				} else {
					Orientation[] ori = parseOrientations(lsplit[12]);
					int ebases = Integer.parseInt(lsplit[8]);
					return new JRecord(
							lsplit[2].replace("chr", ""), 
							Long.parseLong(lsplit[3]), 
							lsplit[5].replace("chr", ""), 
							Long.parseLong(lsplit[6]), 
							ori[0],
							ori[1],
							ebases, caller);
				}
			}
		} catch (IOException e) {
			System.err.println("Error reading file, terminating read, please check integrity of data.");
			e.printStackTrace();
			return null;
		}
	}
	
	@Override
	public boolean hasNext() {
		return next != null;
	}

	@Override
	public JRecord next() {
		JRecord toreturn = next;
		next = parseLine();
		return toreturn;
	}
	
	@Override
	public void close(){
		if (reader != null)
			try {reader.close();} catch (IOException ignore) {}
	}

}
