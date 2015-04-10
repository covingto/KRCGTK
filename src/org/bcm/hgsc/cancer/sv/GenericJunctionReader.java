package org.bcm.hgsc.cancer.sv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.bcm.hgsc.cancer.utils.Orientation;

public class GenericJunctionReader implements JunctionSourceReader {

	private final BufferedReader reader;
	private JRecord next;
	private final Caller caller;
	
	public GenericJunctionReader(String source, String params, int buffer, Caller caller) throws FileNotFoundException{
		this.caller = caller;
		this.reader = new BufferedReader(new FileReader(new File(source)));
		this.next = parseLine();
	}
	
	public GenericJunctionReader(String source, Caller caller) throws FileNotFoundException{
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
			} else if (line.startsWith("@")){
				return parseLine();
			} else {
				String[] lsplit = line.split("\t");
				if (lsplit.length < 8){
					System.err.println("Line did not have 8 fields");
					System.err.print(line);
					return parseLine();
				} else {
					final int bases = Integer.parseInt(lsplit[6]);
					if ((bases < AppendSVGraph.minSVSize) && (lsplit[0].equals(lsplit[2]))){ 
						// not an indel by our standards
						return parseLine();
					}
					return new JRecord(
							lsplit[0], 
							Long.parseLong(lsplit[1]), 
							lsplit[2], 
							Long.parseLong(lsplit[3]), 
							Orientation.valueOf(lsplit[4]),
							Orientation.valueOf(lsplit[5]),
							lsplit[7].equals("DEL") ? 0 : bases, caller);
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
