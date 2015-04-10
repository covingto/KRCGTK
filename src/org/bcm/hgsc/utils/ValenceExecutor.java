package org.bcm.hgsc.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class ValenceExecutor {
	// ===========
	public static String tmpDir = System.getProperty("user.dir");
	public static String scriptDir;
	public static OutputWriter writer;
	public static ExecutorService pool;
	public static List<ValenceAction> active = new ArrayList<ValenceAction>();
	public static Set<String> complete = new HashSet<String>();
	// ===========

	class OutputWriter{
		private BufferedWriter writer;

		OutputWriter(File output) throws IOException{
			this.writer = new BufferedWriter(new FileWriter(output));
		}

		OutputWriter(){	}
		
		public synchronized void writeStreamGobbler(StreamGobbler out, StreamGobbler err, String header){
			write(header);
			write("STD out;");
			for (final String line : out.outputs){
				write(line);
			}
			write("STD err;");
			for (final String line : err.outputs){
				write(line);
			}
		}

		public synchronized void write(String line){
			if (writer != null){
				try {
					writer.write(line);
					writer.newLine();
				} catch (IOException e) {
					e.printStackTrace();
				}
			} else {
				System.out.println(line);
			}
		}
		
		public void close(){
			if (writer != null){
				try{
					writer.close();
				} catch (IOException e){
					e.printStackTrace();
				}
			}
		}
	}

	

	public static void main(String[] args){
		Parser parser = new BasicParser();
		Options options = new Options();
		options.addOption("tmpdir", true, "Prints a summary of the graph to stdout");
		options.addOption("threads", true, "How many threads should be executed at once? [8]");
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
			ValenceExecutor.tmpDir = line.getOptionValue("tmpdir", System.getProperty("user.dir"));
			int threads = Integer.parseInt(line.getOptionValue("threads", "8"));
			pool = Executors.newFixedThreadPool(threads);
			System.out.println("Reading XML file and building dependency tree;");
			
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder builder = factory.newDocumentBuilder();
			Document document = builder.parse(xmlFile);
			document.normalize();
			
			// add the Valence steps to the dependency tree
			addProcedure(document.getDocumentElement(), new HashSet<String>(), new HashMap<String, String>());
			poll();
		} catch (ParseException e) {
			e.printStackTrace();
			System.exit(10);
		} catch (SAXException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void addProcedure(Element procedureElement, Set<String> parentDependencies, Map<String, String> parentReplacements){
		// the procedureElement is the procedure element
		
		// set up the variables
		Set<String> dependencies = new HashSet<String>();
		Map<String, String> replacements = new HashMap<String, String>();
		dependencies.addAll(parentDependencies);
		replacements.putAll(parentReplacements);
		
		// populate the variables for this element
		Element setGlobals = getFirstElement(procedureElement, "setGlobals");
		if (setGlobals != null){
			NodeList sgnl = setGlobals.getChildNodes();
			for (int sg = 0; sg < sgnl.getLength(); sg++){
				Node node = sgnl.item(sg);
				if (node.getNodeType() == Node.ELEMENT_NODE){
					replacements.put("${" + node.getNodeName() + "}", node.getTextContent());
				}
			}
		}
		Element depend = getFirstElement(procedureElement, "dependencies");
		if (depend != null){
			NodeList dnl = depend.getElementsByTagName("depend");
			for (int di = 0; di < dnl.getLength(); di++){
				Node node = dnl.item(di);
				if (node.getNodeType() == Node.ELEMENT_NODE){
					dependencies.add(node.getTextContent());
				}
			}
		}
		
		// add the actions
		NodeList steps = procedureElement.getElementsByTagName("step");
		for (int si = 0; si < steps.getLength(); si++){
			Node node = steps.item(si);
			if (node.getNodeType() == Node.ELEMENT_NODE){
				Element step = (Element) node;
				if (step.getAttribute("type").equals("action")){
					try {
						active.add(new ValenceAction(step, step.getAttribute("name"), replacements, dependencies));
					} catch (DOMException e) {
						System.out.println("Failed to generate action for " + step.getAttribute("name"));
						e.printStackTrace();
					} catch (IOException e) {
						System.out.println("Failed to generate action for " + step.getAttribute("name"));
						e.printStackTrace();
					}
				} else if (step.getAttribute("type").equals("procedure")){
					addProcedure(step, dependencies, replacements);
				}
			}
		}
	}
	
	public static Element getFirstElement(Element parent, String name){
		NodeList nl = parent.getElementsByTagName(name);
		for (int i = 0; i < nl.getLength(); i++){
			Node node = nl.item(i);
			if (node.getNodeType() == Node.ELEMENT_NODE){
				return (Element) node;
			}
		}
		return null;
	}
	
	public static synchronized void poll(){
		List<ValenceAction> toActivated = new ArrayList<ValenceAction>();
		for (ValenceAction action : active){
			if (action.dependenciesCleared(complete)){
				toActivated.add(action);
			}
		}
		for (ValenceAction action : toActivated){
			active.remove(action);
			pool.execute(action);
		}
		if (active.size() == 0){
			// all elements are executing
			pool.shutdown();
		}
	}
}
