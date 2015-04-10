package org.bcm.hgsc.cancer.sv;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import org.json.simple.JSONValue;
import org.xml.sax.SAXException;

public class AppendSVGraph {
	static final int minSVSize = 50;
	/**
	 * Reads in an optional svgraph.ser (or creates a new one) and appends new calls to the graph simultaneously adding 
	 * matchups.
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args){
		if (args.length < 2){
			System.err.println("Not enough arguments supplied");
			System.err.println("Useage; graphObject.ser configurationFile.json");
			System.exit(10);
		}
		//File serFile = new File(args[0]);
		File xmlFile = new File(args[0]);
		SVGraph theGraph;
		if (xmlFile.exists()){
			//System.out.println("Reading serialized file.");
			System.out.println("Reading xml file.");
			try {
				//FileInputStream fileIn = new FileInputStream(serFile);
				//ObjectInputStream in = new ObjectInputStream(fileIn);
				//theGraph = (SVGraph) in.readObject();
				theGraph = SVGraph.readXML(xmlFile);
				//in.close();
				//fileIn.close();
				//System.out.println("Read graph with the following callers;");
				theGraph.printCallers();
			} catch(IOException i) {
				i.printStackTrace();
				return;
			} catch (SAXException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return;
			} catch (ParserConfigurationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return;
			}
		} else {
			System.out.println("Generating new serialized file");
			theGraph = new SVGraph();
		}
		File conffile = new File(args[1]);
		try {
			FileReader fr;
			fr = new FileReader(conffile);

			@SuppressWarnings("unchecked")
			Map<String, Object> conf = (Map<String, Object>) JSONValue.parse(fr);
			fr.close();
			buildGraph(theGraph, conf);

			//FileOutputStream fileOut =
			//		new FileOutputStream(serFile);
			//ObjectOutputStream out = new ObjectOutputStream(fileOut);
			//out.writeObject(theGraph);
			//out.close();
			//fileOut.close();
			//System.out.println("Wrote serialization successfully");
			theGraph.writeXML(xmlFile);
			QuerySVGraph.hiveGraph(theGraph, new File("hiveGraph.svg"), null, null);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (TransformerException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unchecked")
	private static void buildGraph(SVGraph theGraph, Map<String, Object> conf) throws Exception {
		List<Object> settings = (List<Object>) conf.get("settings");
		List<String> chrs = (List<String>) conf.get("chrs");
		for (Object s : settings){
			Map<String, Object> m = (Map<String, Object>) s;
			appendData(theGraph, (String) m.get("source"), (String) m.get("caller"), (String) m.get("conf"), (Long) m.get("buffer"), (String) m.get("platform"), chrs);
		}
		System.out.println("Building links");
		System.out.println("Clearing old matches");
		theGraph.clearMatches();
		System.out.println("Beginning matchdown");
		theGraph.matchDown();
		System.err.println("Beginning matchup");
		theGraph.matchUp();

	}

	private static void appendData(SVGraph theGraph, String source,
			String caller, String conf, Long buffer, String platform, List<String> chrs) throws IOException {
		// create tables if they do not exist
		System.out.println("Parsing data;");
		System.out.println("Source: " + source);
		System.out.println("Caller: " + caller);
		System.out.println("Conf: " + conf);
		System.out.println("Platform: " + platform);
		Caller newCaller = new Caller(caller, conf, buffer.intValue(), platform);
		theGraph.addCaller(newCaller);

		// add the caller info to the database
		try {
			JunctionSourceReader jsr;
			if (caller.equals("BreakDown")){
				jsr = new BreakDownJunctionReader(source, newCaller);
			} else if (caller.equals("BreakDancer")){
				jsr = new BreakDancerReader(source, newCaller);
			} else if (caller.equals("PInDel")){
				jsr = new PInDelJunctionReader(source, newCaller);
			} else if (caller.equals("GATK")) {
				jsr = new VCFJunctionReader(source, newCaller);
			} else if (caller.equals("HoneyTails")) {
				jsr = new HoneyTailsJunctionReader(source, newCaller);
			} else if (caller.equals("HoneySpots")) {
				jsr = new HoneySpotsJunctionReader(source, newCaller);
			} else {
				jsr = new GenericJunctionReader(source, newCaller);
			}
			int nrecs = 0;
			while (jsr.hasNext()){
				final JRecord jr = jsr.next();
				if (chrs != null && (!chrs.contains(jr.chra()) || !chrs.contains(jr.chrb())) ){ 
					// System.out.println("Excluede chrs: " + jr.chra() + ", " + jr.chrb());
					continue; }
				theGraph.addJunction(new GraphJunction(jr));
				nrecs++;
			}
			jsr.close();
			System.out.println("Wrote " + nrecs + " to database");
		} catch (IOException readExcept){
			System.err.println("Error reading source file: " + source);
			readExcept.printStackTrace();
			throw readExcept;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
