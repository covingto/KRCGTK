package org.bcm.hgsc.cancer.sv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.bcm.hgsc.cancer.utils.Orientation;

public class VCFJunctionReader implements JunctionSourceReader {

	private final BufferedReader reader;
	private JRecord next;
	private final Caller caller;
	
	public VCFJunctionReader(String source, String params, int buffer, Caller caller) throws FileNotFoundException{
		this.caller = caller;
		this.reader = new BufferedReader(new FileReader(new File(source)));
		this.next = parseLine();
	}
	
	public VCFJunctionReader(String source, Caller caller) throws FileNotFoundException{
		this(source, "default", 1000, caller);
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
				if (lsplit.length < 7){
					System.err.println("Line did not have 7 fields");
					System.err.print(line);
					return parseLine();
				} else {
					//final int ebases = lsplit[5].replace("\"", "").length();
					int dbases = lsplit[3].length();
					int ibases = lsplit[4].length();
					int ebases;
					if (Math.max(lsplit[3].length(), lsplit[4].length()) < AppendSVGraph.minSVSize) { 
						// recurse, bases are not enough to warrant an SV (in our sense)
						return parseLine(); 
					} else if (dbases > ibases) {
						ebases = 0;
					} else {
						ebases = ibases;
					}
					
					return new JRecord(
							lsplit[0], 
							Long.parseLong(lsplit[1]), 
							lsplit[0], 
							Long.parseLong(lsplit[1]) + dbases, 
							Orientation.POS,
							Orientation.POS,
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
