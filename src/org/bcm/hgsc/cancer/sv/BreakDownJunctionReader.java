package org.bcm.hgsc.cancer.sv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.bcm.hgsc.cancer.utils.Orientation;

public class BreakDownJunctionReader implements JunctionSourceReader {
	private final BufferedReader reader;
	private JRecord next;
	private final Caller caller;
	
	public BreakDownJunctionReader(String source, String params, int buffer, Caller caller) throws FileNotFoundException{
		this.caller = caller;
		this.reader = new BufferedReader(new FileReader(new File(source)));
		this.next = parseLine();
	}
	
	public BreakDownJunctionReader(String source, Caller caller) throws FileNotFoundException{
		this(source, "default", 1000, caller);
	}
	
	private static Orientation[] parseOrientation(String s){
		boolean pp = s.contains("POS,POS");
		boolean pn = s.contains("POS,NEG");
		boolean np = s.contains("NEG,POS");
		if (pp && !(pn || np)){
			// only positive, easy case
			return new Orientation[] {Orientation.POS, Orientation.POS};
		} else if (pn && np){
			// most complicated case, reciprocal translocation
			return new Orientation[] {Orientation.RECIP, Orientation.RECIP};
		} else if (pp && (pn || np)){
			// can't really tell what is happening here, two kinds of orientaitons supported
			return new Orientation[] {Orientation.AMBIGUOUS, Orientation.AMBIGUOUS};
		} else if (pn){
			return new Orientation[] {Orientation.POS, Orientation.NEG};
		} else if (np){
			return new Orientation[] {Orientation.NEG, Orientation.POS};
		} else {
			return new Orientation[] {Orientation.AMBIGUOUS, Orientation.AMBIGUOUS};
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
					Orientation[] ori = parseOrientation(lsplit[8]);
					
					int rebases = Integer.parseInt(lsplit[7]);
					int ebases = rebases > 0 ? rebases : 0;
					return new JRecord(
							lsplit[1].replace("chr", ""), 
							Long.parseLong(lsplit[2]), 
							lsplit[3].replace("chr", ""), 
							Long.parseLong(lsplit[4]), 
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
