package org.bcm.hgsc.cancer.sv;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
import org.w3c.dom.Document;
import org.w3c.dom.DocumentType;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

public class QuerySVGraph {
	// ==================================
	public static class SVGGraphJunctionHolder{
		private final GraphJunction j;
		private final double x;
		private final double y;
		private final double theta;
		private final int callerInt;
		
		public SVGGraphJunctionHolder(double x, double y, double theta, int callerInt, GraphJunction j){
			this.j = j;
			this.x = x;
			this.y = y;
			this.theta = theta;
			this.callerInt = callerInt;
		}
		
		public double getX(){
			return this.x;
		}
		
		public int getCallerInt(){
			return this.callerInt;
		}
		
		public double getY(){
			return this.y;
		}
		
		public GraphJunction getJ(){
			return this.j;
		}
		
		public double getTheta(){
			return this.theta;
		}
	}
	
	// ==================================
	public enum HiveType{
		ByCaller, ByClass, ByChr;
	}
	
	// ==================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("summary", false, "Prints a summary of the graph to stdout");
		options.addOption("svghive", true, "Generates a hive plot of the graph.  Must supply a configuration file.");
		options.addOption("familyCount", false, "Generates a table of cluster sizes from the graph.");
		options.addOption("callerSummary", false, "Generates a summary of callers from the graph.");
		
		// filtering and subsetting options
		options.addOption("filterHit", true, "Filters hit junctions for the given caller");
		options.addOption("filterNoHit", true, "Filters non-hit junctions for the given caller, mutually exclusive with filterHit");
		options.addOption("fullHit", false, "Flag for filterHit indicating if full or partial hits should be considered");
		options.addOption("filterMask", true, "List of callers to mask when considering the filter for filterHit");
		options.addOption("filterOutput", true, "File to write the filter results (Generic Junction format) to [stdout]");
		options.addOption("writeJunctions", false, "Write junctions out to a file");
		// options.addOption()
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line;
		String[] reqargs;
		File xmlFile = null;
		try {
			line = parser.parse(options, args);
			reqargs = line.getArgs();
			if (reqargs.length < 1){ throw new ParseException("Required arguments not found!"); }
			xmlFile = new File(reqargs[0]);
			if (!xmlFile.exists()){ throw new ParseException("XML file not found!"); }

			SVGraph theGraph;
			System.out.println("Reading serialized file.");

			//FileInputStream fileIn = new FileInputStream(serFile);
			//ObjectInputStream in = new ObjectInputStream(fileIn);
			theGraph = SVGraph.readXML(xmlFile);
			//in.close();
			//fileIn.close();
			//theGraph.update();
			System.out.println("Read graph with the following callers;");
			theGraph.printCallers();
			if (line.hasOption("familyCount")){
				QuerySVGraph.printFamilyCount(theGraph, new File("familyCounts.txt"));
			}
			if (line.hasOption("summary")){
				QuerySVGraph.summary(theGraph, System.out);
			}
			if (line.hasOption("svghive")){
				QuerySVGraph.hiveGraph(theGraph, new File("hivePlot.svg"), HiveType.ByCaller, null);
			}
			if (line.hasOption("callerSummary")){
				QuerySVGraph.printCallerSummary(theGraph);
			} 
			if (line.hasOption("filterHit")){
				Caller filterCaller = theGraph.getCallerByName(line.getOptionValue("filterHit"));
				Set<String> mask = new HashSet<String>();
				if (line.hasOption("filterMask")){
					String maskList = line.getOptionValue("filterMask");
					String[] maskL = maskList.split(",");
					Map<String,String> callerLookup = new HashMap<String, String>();
					for (Caller c : SVGraph.callers.values()){
						callerLookup.put(c.getName(), c.getUUID());
					}
					for (int m = 0; m < maskL.length; m++){
						String toAdd = callerLookup.get(maskL[m]);
						if (toAdd != null){
							mask.add(toAdd);
							System.err.println("Adding " + maskL[m] + " to mask, uuid=" + toAdd);
						} else {
							System.err.println("Can't find caller " + maskL[m]);
						}
					}
				}
				if (filterCaller == null){ 
					System.err.println("Could not find caller " + line.getOptionValue("filterHit") + " in caller list!");
				} else {
					printHitJunctions(filterCaller.getUUID(), line.hasOption("fullHit"), mask, line.getOptionValue("filterOutput", null));
				}
			} else if (line.hasOption("filterNoHit")){
				Caller filterCaller = theGraph.getCallerByName(line.getOptionValue("filterNoHit"));
				Set<String> mask = new HashSet<String>();
				if (line.hasOption("filterMask")){
					String maskList = line.getOptionValue("filterMask");
					String[] maskL = maskList.split(",");
					Map<String,String> callerLookup = new HashMap<String, String>();
					for (Caller c : SVGraph.callers.values()){
						callerLookup.put(c.getName(), c.getUUID());
					}
					for (int m = 0; m < maskL.length; m++){
						String toAdd = callerLookup.get(maskL[m]);
						if (toAdd != null){
							mask.add(toAdd);
							System.err.println("Adding " + maskL[m] + " to mask, uuid=" + toAdd);
						} else {
							System.err.println("Can't find caller " + maskL[m]);
						}
					}
				}
				if (filterCaller == null){ 
					System.err.println("Could not find caller " + line.getOptionValue("filterNoHit") + " in caller list!");
				} else {
					printHitJunctions(filterCaller.getUUID(), line.hasOption("fullHit"), mask, line.getOptionValue("filterOutput", null));
				}
			} else if (line.hasOption("writeJunctions")){
				Set<String> mask = new HashSet<String>();
				if (line.hasOption("filterMask")){
					String maskList = line.getOptionValue("filterMask");
					String[] maskL = maskList.split(",");
					Map<String,String> callerLookup = new HashMap<String, String>();
					for (Caller c : SVGraph.callers.values()){
						callerLookup.put(c.getName(), c.getUUID());
					}
					for (int m = 0; m < maskL.length; m++){
						String toAdd = callerLookup.get(maskL[m]);
						if (toAdd != null){
							mask.add(toAdd);
							System.err.println("Adding " + maskL[m] + " to mask, uuid=" + toAdd);
						} else {
							System.err.println("Can't find caller " + maskL[m]);
						}
					}
				} 
				printJunctions(mask, line.getOptionValue("filterOutput", null));
			}
		} catch (ParseException e) {
			e.printStackTrace();
			formatter.printHelp("java -jar QuerySVGraph.jar [options] [serialized graph]", options);
			System.exit(10);
		} catch(IOException i) {
			i.printStackTrace();
			System.exit(10);
		}catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (TransformerException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SAXException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// ==================================
	public static void printFamilyCount(SVGraph g, File outfile) throws IOException{
		Integer[] fc = g.generateJunctionFamilies();
		BufferedWriter writer = new BufferedWriter(new FileWriter(outfile));
		writer.write("Index\tFamilyCount");
		writer.newLine();
		for (int i = 0; i < fc.length; i++){
			writer.write(i + "\t" + fc[i]);
			writer.newLine();
		}
		writer.close();
	}
	
	// ==================================
	public static void printCallerSummary(SVGraph g){
		Map<String, List<GraphJunction>> lJunctionsByCaller = g.getJunctionsByCaller();
		System.out.print("Caller\tJunctions\tInterchromosomal\tIntrachromosomal\tCoCalled");
		for (final String cid : lJunctionsByCaller.keySet()){
			final Caller caller = g.getCaller(cid);
			System.out.print("\t"+caller.getName());
		}
		System.out.print("\n");
		for (final String cid : lJunctionsByCaller.keySet()){
			final Caller caller = g.getCaller(cid);
			final String callerName = caller.getName();
			// loop through the junctions within the caller and do some counting
			List<GraphJunction> calls = lJunctionsByCaller.get(cid);
			final int nCalls = calls.size();
			// loop through the nodes and add the callers.
			Map<String, List<GraphJunction>> corespondingJunctions = new HashMap<String, List<GraphJunction>>();
			for (String ocid : lJunctionsByCaller.keySet()){
				if (ocid == cid){ continue; }
				corespondingJunctions.put(ocid, new ArrayList<GraphJunction>());
			}
			int coCalled = 0;
			int interChromCount = 0;
			int intraChromCount = 0;
			for (GraphJunction j : calls){
				for (String jc : j.partnerCallers){
					corespondingJunctions.get(jc).add(j);
				}
				if (j.chra().equals(j.chrb())){
					intraChromCount++;
				} else {
					interChromCount++;
				}
				if (j.partnerCallers.size() > 1){
					coCalled++;
				}
			}
			String outstr = callerName + "\t" + nCalls + "\t" + interChromCount + "\t" + intraChromCount + "\t" + coCalled;
			for (final String ocid : lJunctionsByCaller.keySet()){
				if (ocid == cid) { outstr = outstr + "\tNA"; }
				else {
					outstr = outstr + "\t" + corespondingJunctions.get(ocid).size();
				}
			}
			System.out.println(outstr);
		}
	}
	
	public static boolean hasFullHit(GraphJunction j, Set<String> mask){
		for (String pj : j.partners){
			GraphJunction p = SVGraph.junctions.get(pj);
			if (mask.contains(p.getCallerUUID())){ continue; }
			switch (j.hitType(p)){
			case Full:
				return true;
			case None:
				break;
			case Right:
			case Left:
				break;
			default:
				break;
			}
		}
		return false;
	}
	
	public static boolean hasAnyHit(GraphJunction j, Set<String> mask){
		for (String pj : j.partners){
			GraphJunction p = SVGraph.junctions.get(pj);
			if (mask.contains(p.getCallerUUID())){ continue; }
			return true;
		}
		return false;
	}
	
	public static void printNonHitJunctions(String caller, boolean fullHit, Set<String> mask, String output) throws IOException{
		BufferedWriter writer;
		if (output != null){
			writer = new BufferedWriter(new PrintWriter(new FileWriter(output, true)));
		} else {
			writer = new BufferedWriter(new OutputStreamWriter(System.out));
		}
		for (GraphJunction j : SVGraph.junctions.values()){
			if (j.getCallerUUID().equals(caller)){
				if (fullHit && !hasFullHit(j, mask)){
					j.printJunction(writer); // here we have to cast the null as a string to avoid a conflict
				} else if (!hasAnyHit(j, mask)){
					j.printJunction(writer);
				}
			}
		}
		writer.close();
	}
	
	public static void printJunctions(Set<String> mask, String output) throws IOException{
		BufferedWriter writer;
		if (output != null){
			writer = new BufferedWriter(new PrintWriter(new FileWriter(output, true)));
		} else {
			writer = new BufferedWriter(new OutputStreamWriter(System.out));
		}
		for (GraphJunction j : SVGraph.junctions.values()){
			if (mask.contains(j.getCallerUUID())){ continue; }
			else { j.printJunction(writer); }
		}
		writer.close();
	}
	
	public static void printHitJunctions(String caller, boolean fullHit, Set<String> mask, String output) throws IOException{
		BufferedWriter writer;
		if (output != null){
			writer = new BufferedWriter(new PrintWriter(new FileWriter(output, true)));
		} else {
			writer = new BufferedWriter(new OutputStreamWriter(System.out));
		}
		for (GraphJunction j : SVGraph.junctions.values()){
			if (j.getCallerUUID().equals(caller)){
				if (fullHit && hasFullHit(j, mask)){
					j.printJunction(writer); // here we have to cast the null as a string to avoid a conflict
				} else if (hasAnyHit(j, mask)){
					j.printJunction(writer);
				}
			}
		}
		writer.close();
	}

	// ==================================
	/**
	 * Generates a hive plot of the graph and writes to the specified file.  Additional options can be passed using the config file.
	 * @param g
	 * @param outfile
	 * @param type
	 * @param config
	 * @throws ParserConfigurationException 
	 * @throws TransformerException 
	 */
	public static void hiveGraph(SVGraph g, File outfile, HiveType type, File config) throws ParserConfigurationException, TransformerException{
		System.out.println("Generating graph subsets");
		final double pointR = 200.0;
		final double outerConR = pointR * 1.25;
		final int numJunctions = g.getJunctions().size();
		
		Map<String, List<GraphJunction>> lJunctionsByCaller = g.getJunctionsByCaller();
		Map<String, SVGGraphJunctionHolder> lJunctionHolderMap = new HashMap<String, SVGGraphJunctionHolder>();

		int ncallers = lJunctionsByCaller.size(); // we need to separate the coordinate system into sections of this size
		double degreesPerBuffer	= ((2 * Math.PI) / 360) * 6; 
		double degreesPerObs	= ((2 * Math.PI) - (degreesPerBuffer * ncallers)) / numJunctions; 
		//double degreesPerCaller = (2 * Math.PI) / ncallers;
		//double degreesPerArc 	= degreesPerCaller * 0.9;
		//double degreesPerBuffer = degreesPerCaller * 0.05;
		
		// lJunctionsByCaller are now organized in such a way that the junctions are clustered by confirmation number
		// we can now do some work
		
		// set up the document
		
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		DocumentType docType = builder.getDOMImplementation().createDocumentType("svg", "-//W3C//DTD SVG 1.1//EN", "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd");
		Document document = builder.getDOMImplementation().createDocument("http://www.w3.org/2000/svg", "svg", docType);
		document.normalize();
		Element svg = document.getDocumentElement();
		svg.setAttribute("width", "1600px");
		svg.setAttribute("height", "1600px");
		svg.setAttribute("viewBox", "0 0 1600 1600");
		// make the style element
		Element defs = document.createElement("defs");
		Element style = document.createElement("style");
		defs.appendChild(style);
		svg.appendChild(defs);
		style.setAttribute("type", "text/css");
		// set up the color transitions
		String callerStyle = "";
		for (int i = 0; i < ncallers; i++){
			final double p = (double) i / (double) ncallers;
			final Double R = Color.YELLOW.getRed() * p + Color.BLUE.getRed() * ( 1 - p );
			final Double G = Color.YELLOW.getGreen() * p + Color.BLUE.getGreen() * ( 1 - p );
			final Double B = Color.YELLOW.getBlue() * p + Color.BLUE.getBlue() * ( 1 - p );
			callerStyle += "circle.hitCount" + i + " {" +
					"	fill: rgb(" + R.intValue() + "," + G.intValue() + "," + B.intValue() + "); }\n";
		}
		style.appendChild(document.createCDATASection(callerStyle +
				"path { " +
				"	fill-opacity:0; stroke-opacity: 0.1; }\n" +
				"path.Partial { " +
				"	stroke: red; }\n" +
				"path.Full { " +
				"	stroke: blue; }\n"));
		
		Element graphArea = document.createElement("g");
		graphArea.setAttribute("transform", "translate(800 800)");
		Element graphNodesGroup = document.createElement("g");
		graphNodesGroup.setAttribute("class", "grapharea");
		graphNodesGroup.setAttribute("id", "nodesGroup");
		
		Element graphEdgesGroup = document.createElement("g");
		graphEdgesGroup.setAttribute("class", "grapharea");
		graphEdgesGroup.setAttribute("id", "edgesGroup");
		
		Element titleGroup = document.createElement("g");
		titleGroup.setAttribute("class", "titlearea");
		titleGroup.setAttribute("id", "titleGroup");
		
		graphArea.appendChild(graphEdgesGroup);
		graphArea.appendChild(graphNodesGroup);
		graphArea.appendChild(titleGroup);
		
		svg.appendChild(graphArea);
		// the graph consists of lines representing the callers
		
		// loop across the junctions and start to build the hive
		int currentCallerIndex = 0;
		double rollingOffset = degreesPerBuffer / 2; // the rolling offset will be incremented
		for (String c : lJunctionsByCaller.keySet()){
			final List<GraphJunction> flJunctions = lJunctionsByCaller.get(c); // these are in sorted order
			final int lsize = flJunctions.size();
			System.out.println("Generating data for caller " + g.getCaller(c).getName());
			
			
			int jindex = 0;
			for (GraphJunction j : flJunctions){
				final Set<String> partners = j.partners;
				final Caller jCaller = g.getCaller(j.getCallerUUID());
				final double theta =  (jindex * degreesPerObs) + rollingOffset;
				final double x = polarToRectangularX(pointR, theta);
				final double y = polarToRectangularY(pointR, theta);
				// add the node to the graph
				final Element thisNode = document.createElement("circle");
				thisNode.setAttribute("cx", String.valueOf(x));
				thisNode.setAttribute("cy", String.valueOf(y));
				thisNode.setAttribute("r", "5");
				thisNode.setAttribute("class", "node " + jCaller.getName() + " hitCount" + partners.size());
				thisNode.setAttribute("id", j.getUUID());
				
				graphNodesGroup.appendChild(thisNode);
				for (String pKey : partners){
					SVGGraphJunctionHolder n = lJunctionHolderMap.get(pKey);
					if (n != null){
						final Element edge = document.createElement("path");
						final Caller pCaller = g.getCaller(n.j.getCallerUUID());
						switch (j.hitType(n.j)){
						case Left:
						case Right:
							// figure out if the arc is forward or backward
							// arc 1 will have been drawn before arc 2, if the difference between theta1 and theta2 (this theta) is < pi then use positive angles, else negative
							// org; String.valueOf(currentCallerIndex - n.callerInt < (ncallers / 2) ? 0 : 1)
							edge.setAttribute("d", 
									"M " + String.valueOf(x) + " " + String.valueOf(y) + 
									" A " + outerConR + " " + outerConR + " 0 1 " + String.valueOf(theta - n.theta < Math.PI ? 0 : 1) + 
									String.valueOf(n.x) + " " + String.valueOf(n.y));
							edge.setAttribute("class", "Partial C" + j.getCallerUUID() + " C" + n.j.getCallerUUID() + " " + jCaller.getPlatform() + " " + pCaller.getPlatform());
							break;
						case Full:
							edge.setAttribute("d", 
									"M " + String.valueOf(x) + " " + String.valueOf(y) +
									" Q 0 0 " + String.valueOf(n.x) + " " + String.valueOf(n.y));
							edge.setAttribute("class", "Full C" + j.getCallerUUID() + " C" + n.j.getCallerUUID() + " " + jCaller.getPlatform() + " " + pCaller.getPlatform());
							break;
						default:
							continue;
						}
						graphEdgesGroup.appendChild(edge);
					}
				}
				// add the node to the lJunctionHolderMap
				lJunctionHolderMap.put(j.getUUID(), new SVGGraphJunctionHolder(x, y, theta, currentCallerIndex, j));
				jindex++;
			}
			currentCallerIndex++;
			double newOffset = (degreesPerObs * lsize) + degreesPerBuffer + rollingOffset;
			final double thetaMiddle = (rollingOffset + newOffset) / 2;
			final double xLabel = polarToRectangularX(pointR * 2, thetaMiddle);
			final double yLabel = polarToRectangularY(pointR * 2, thetaMiddle);
			final Element callerLab = document.createElement("text");
			callerLab.setAttribute("x", String.valueOf(xLabel));
			callerLab.setAttribute("y", String.valueOf(yLabel));
			titleGroup.appendChild(callerLab);
			final Caller thisCaller = g.getCaller(c);
			callerLab.appendChild(document.createTextNode("Caller " + (currentCallerIndex + 1) + ": " + thisCaller.getName()));
			System.out.println("moving from " + rollingOffset + " to " + newOffset);
			rollingOffset = newOffset;
//			Element lastPos = document.createElement("text");
//			lastPos.setAttribute("x", String.valueOf(lastX));
//			lastPos.setAttribute("y", String.valueOf(lastY));
//			titleGroup.appendChild(callerLab);
//			lastPos.appendChild(document.createTextNode(String.valueOf(currentCallerIndex)));
			
		}
		
		// stream out the output
		TransformerFactory tFactory = TransformerFactory.newInstance();
		Transformer transformer = tFactory.newTransformer();
		transformer.setOutputProperty(OutputKeys.INDENT, "yes");
		transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
		DOMSource source = new DOMSource(document);
		StreamResult result = new StreamResult(outfile);
		transformer.transform(source, result);
	}
	
	// ==================================
//	private static double scaleThetaPlacement(int lrank, int lsize, double degreesPerArc){
//		return degreesPerArc * (double) lrank / (double) lsize;
//	}
	
	// ==================================
	private static double polarToRectangularX(double r, double t){
		return r * Math.cos(t);
	}
	
	// ==================================
	private static double polarToRectangularY(double r, double t){
		return r * Math.sin(t);
	}

	// ==================================
	public static void summary(SVGraph g, OutputStream os) throws IOException{
		// print the number of nodes
		os.write(("Nodes: " + g.nodeCount() + String.format("%n")).getBytes());
		os.write(("Edges: " + g.edgeCount() + String.format("%n")).getBytes());
		os.write(String.format("%n").getBytes());
		//build caller information
		// map of uuid of the caller (for lookup) and a list of junctions for that caller
		Map<String, List<GraphJunction>> callerhits = g.getCallerData();
		for (String k : callerhits.keySet()){
			final List<GraphJunction> hits = callerhits.get(k);
			Map<String, Integer> partners = new HashMap<String, Integer>();
			for (GraphJunction gj : hits){
				Set<String> partnerCallers = g.getPartnerCallers(gj);
				for (String pc : partnerCallers){
					if (!partners.containsKey(pc)){ partners.put(pc, 0); }
					partners.put(pc, partners.get(pc) + 1);
				}
			}
			os.write(("[" + g.getCaller(k).getName() + ", " + k + "]\t" + hits.size() ).getBytes());
			for (String p : partners.keySet()){
				os.write(("\t(" + g.getCaller(p).getName() + ", " + partners.get(p) + ")").getBytes());
			}
			os.write(String.format("%n").getBytes());
		}
		// print; [callername, caller int] \t number of calls \t (cocaller 1, hits) \t (cocaller 2, hits) ...

	}

}
