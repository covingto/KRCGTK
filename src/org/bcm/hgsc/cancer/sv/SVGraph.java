package org.bcm.hgsc.cancer.sv;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.bcm.hgsc.cancer.sv.GraphJunction.HitType;
import org.bcm.hgsc.cancer.utils.Orientation;
import org.bcm.hgsc.utils.Settings;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class SVGraph implements Serializable{

	/**
	 * serial version as of 14 March 2014
	 */
	private static final long serialVersionUID = 1L;
	//private static Map<String, GraphJunction> gjunctions;
	//private static boolean globalset = false;
	public static Map<String, GraphJunction> junctions = new HashMap<String, GraphJunction>();
	private Map<String, List<String>> aSortedJunctions = new HashMap<String, List<String>>();
	private Map<String, List<String>> bSortedJunctions = new HashMap<String, List<String>>();
	public static Map<String, Caller> callers = new HashMap<String, Caller>();
	private boolean sorted = false;
	private static boolean familiesMade = false;
	private static Integer[] familyCounts = null;

	private class MergerWorker implements Runnable {
		private final List<String> gjl;
		MergerWorker(List<String> gjl){ this.gjl = gjl; }
		@Override
        public void run(){
			System.err.println("Matching up " + gjl.size() + " records");
			for (String js : gjl){
				GraphJunction j = junctions.get(js);
				for (GraphJunction p : getMatchingJunctions(j)){
					j.addPartner(p.getUUID());
					j.addPartnerCaller(p.getCallerUUID());
					p.addPartner(j.getUUID());
					p.addPartnerCaller(j.getCallerUUID());
				}
			}
		}
	}

	private class MatchDownWorker implements Runnable{
		private final List<GraphJunction> gjl;
		private final String ck;
		MatchDownWorker(List<GraphJunction> list, String ck){
			this.gjl = list;
			this.ck = ck;
		}
		@Override
        public void run(){
			//System.out.println("Matching down junctions from caller; " + SVGraph.callers.get(ck).getName());
			List<GraphJunction> callerJunctionsASorted = new ArrayList<GraphJunction>();
			List<GraphJunction> callerJunctionsBSorted = new ArrayList<GraphJunction>();
			List<GraphJunction> junctions = this.gjl;
			callerJunctionsASorted.addAll(this.gjl);
			callerJunctionsBSorted.addAll(this.gjl);
			Collections.sort(callerJunctionsASorted, GraphJunction.aComparator);
			Collections.sort(callerJunctionsBSorted, GraphJunction.bComparator);
			Set<GraphJunction> seenSet = new HashSet<GraphJunction>();
			List<GraphJunction> addingJunctions = new ArrayList<GraphJunction>();
			List<GraphJunction> toRemove = new ArrayList<GraphJunction>();
			int removedCount = 0;
			for (GraphJunction j : junctions){
				// find the family for this junction
				if (seenSet.contains(j)){ continue; }
				Set<GraphJunction> newSet = new HashSet<GraphJunction>();
				newSet.add(j);
				_popMatchDownGJSet(newSet, j, callerJunctionsASorted, callerJunctionsBSorted);
				seenSet.addAll(newSet);
				if (newSet.size() > 1){
					toRemove.addAll(newSet); // we need to remove all of this family of junctions
					//System.out.println("Matchdown removing " + newSet.size() + " junctions.");
					removedCount++;
					GraphJunction newJunction = _combineGraphJunctions(newSet);
					addingJunctions.add(newJunction);
				}
			}
			for (GraphJunction rj : toRemove){
				SVGraph.junctions.remove(rj.getUUID());
			}
			for (GraphJunction aj : addingJunctions){
				try {
					addJunction(aj);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			System.out.println("Matching down junctions from caller; " + SVGraph.callers.get(this.ck).getName() + " removed " + removedCount + " junctions");
		}
	}

	public class JunctionMerger{

	}
	// =====================
	public static final Comparator<String> aComparator = new Comparator<String>(){

		@Override
		public int compare(String o1, String o2) {
			GraphJunction o1j = SVGraph.junctions.get(o1);
			GraphJunction o2j = SVGraph.junctions.get(o2);
			if (o1j.chra().equals(o2j.chra())){
				return (int) (o1j.posa() - o2j.posa());
			} else {
				return o1j.chra().compareTo(o2j.chra());
			}
		}

	};
	// =====================
	public static final Comparator<String> bComparator = new Comparator<String>(){

		@Override
		public int compare(String o1, String o2) {
			GraphJunction o1j = SVGraph.junctions.get(o1);
			GraphJunction o2j = SVGraph.junctions.get(o2);
			if (o1j.chrb().equals(o2j.chrb())){
				return (int) (o1j.posb() - o2j.posb());
			} else {
				return o1j.chrb().compareTo(o2j.chrb());
			}
		}

	};
	// =====================
	public void addJunction(GraphJunction j) throws Exception{
		if (!callers.containsKey(j.getCallerUUID())){
			throw new Exception("Caller not in callers list.");
		}
		junctions.put(j.getUUID(), j);
		sorted = false;
	}


	// =====================
	public Map<String, List<GraphJunction>> getJunctionsByCaller(){
		Map<String, List<GraphJunction>> lJunctionsByCaller = new HashMap<String, List<GraphJunction>>();
		for (GraphJunction j : getJunctions()){
			final String caller = j.getCallerUUID();
			if (!lJunctionsByCaller.containsKey(caller)){ lJunctionsByCaller.put(caller, new ArrayList<GraphJunction>()); }
			lJunctionsByCaller.get(caller).add(j);
		}
		for (String k : lJunctionsByCaller.keySet()){
			Collections.sort(lJunctionsByCaller.get(k), GraphJunction.interactionCountComparator);
		}
		return lJunctionsByCaller;
	}


	// =====================
	public Integer[] generateJunctionFamilies(){
		return generateJunctionFamilies(false);
	}

	// =====================
	public Integer[] generateJunctionFamilies(boolean override){
		if (!override && familiesMade && familyCounts != null){ return familyCounts; }
		Set<GraphJunction> seen = new HashSet<GraphJunction>(SVGraph.junctions.size());
		Map<String, GraphJunction> junctions = SVGraph.junctions;
		List<Integer> seenCounts = new ArrayList<Integer>();
		int familyInt = 0;
		for (GraphJunction j : junctions.values()){
			if (seen.contains(j)){ continue; }
			Set<GraphJunction> familyMembers = new HashSet<GraphJunction>();
			familyMembers.add(j);
			_addPartners(familyMembers, j);
			seen.addAll(familyMembers);
			seenCounts.add(familyMembers.size());
			for (GraphJunction fj : familyMembers){
				fj.setFamilyInt(familyInt);
				fj.setFamilyCount(familyMembers.size());
			}
		}
		familyCounts = seenCounts.toArray(new Integer[seenCounts.size()]);
		familiesMade = true;
		return familyCounts;
	}

	// =====================
	private static void _addPartners(Set<GraphJunction> familyMembers, GraphJunction j) {
		for (String p : j.partners){
			final GraphJunction pj = SVGraph.junctions.get(p);
			if (pj != null){
				if (familyMembers.contains(pj)){ continue; } // don't add a member that is already in the set
				familyMembers.add(pj);
				// recurse
				_addPartners(familyMembers, pj);
			} else {
				System.err.println("No junction record not found in graph; " + p);
			}
		}

	}

	// =====================
	public void addCaller(Caller c){
		callers.put(c.getUUID(), c);
	}

	// =====================
	public void writeXML(String destination) throws ParserConfigurationException, TransformerException{
		writeXML(new File(destination));
	}

	// ==================================
	public void writeXML(File destination) throws ParserConfigurationException, TransformerException{
		// Generate an XML document and load the data
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document document = builder.newDocument();
		document.normalize();
		Element graph = document.createElement("graph");
		document.appendChild(graph);
		Element callersElement = document.createElement("callers");
		Element nodesElement = document.createElement("nodes");
		graph.appendChild(callersElement);
		graph.appendChild(nodesElement);
		// add the callers
		for (String ck : callers.keySet()){
			Caller thisC = callers.get(ck);
			Element thisCNode = document.createElement("caller");
			thisCNode.setAttribute("uuid", thisC.getUUID());
			thisCNode.setAttribute("name", thisC.getName());
			thisCNode.setAttribute("params", thisC.getParams());
			thisCNode.setAttribute("buffer", String.valueOf(thisC.getBuffer()));
			thisCNode.setAttribute("platform", thisC.getPlatform());
			callersElement.appendChild(thisCNode);
		}
		for (String jk : junctions.keySet()){
			GraphJunction thisJ = junctions.get(jk);
			Element thisJNode = document.createElement("junction");
			thisJNode.setAttribute("uuid", thisJ.getUUID());
			thisJNode.setAttribute("caller", thisJ.getCallerUUID());
			thisJNode.setAttribute("chra", thisJ.chra());
			thisJNode.setAttribute("chrb", thisJ.chrb());
			thisJNode.setAttribute("posa", String.valueOf(thisJ.posa()));
			thisJNode.setAttribute("posb", String.valueOf(thisJ.posb()));
			thisJNode.setAttribute("oapos", String.valueOf(thisJ.oapos()));
			thisJNode.setAttribute("obpos", String.valueOf(thisJ.obpos()));
			thisJNode.setAttribute("ebases", String.valueOf(thisJ.ebases()));
			thisJNode.setAttribute("familyInt", String.valueOf(thisJ.getFamilyInt()));
			thisJNode.setAttribute("familyCount", String.valueOf(thisJ.getFamilyCount()));
			for (String p : thisJ.partners){
				Element partnerNodeElement = document.createElement("partner");
				partnerNodeElement.setAttribute("id", p);
				partnerNodeElement.setAttribute("caller", SVGraph.junctions.get(p).getCallerUUID());
				thisJNode.appendChild(partnerNodeElement);
			}
			nodesElement.appendChild(thisJNode);

		}
		// stream out the output
		TransformerFactory tFactory = TransformerFactory.newInstance();
		Transformer transformer = tFactory.newTransformer();
		transformer.setOutputProperty(OutputKeys.INDENT, "yes");
		transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
		DOMSource source = new DOMSource(document);
		StreamResult result = new StreamResult(destination);
		transformer.transform(source, result);
	}

	// ==================================
	public static SVGraph readXML(String file) throws Exception{
		return readXML(new File(file));
	}

	// ==================================
	public static SVGraph readXML(File file) throws Exception{
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document document = builder.parse(file);
		document.normalize();

		// start building the graph
		SVGraph theGraph = new SVGraph();

		// get the callers
		Element callersElement = (Element) document.getElementsByTagName("callers").item(0);
		NodeList callerNodesList = callersElement.getElementsByTagName("caller");
		for (int i = 0; i < callerNodesList.getLength(); i++){
			Node node = callerNodesList.item(i);
			if (node.getNodeType() == Node.ELEMENT_NODE){
				Element callerElement = (Element) node;
				final String name = callerElement.getAttribute("name");
				final Integer buffer = Integer.parseInt(callerElement.getAttribute("buffer"));
				final String params = callerElement.getAttribute("params");
				final String uuid = callerElement.getAttribute("uuid");
				final String platform = callerElement.getAttribute("platform");
				theGraph.addCaller(new Caller(name, params, buffer, platform, uuid));
			}
		}

		// get the nodes
		Element nodesElement = (Element) document.getElementsByTagName("nodes").item(0);
		NodeList nodesNodeList = nodesElement.getElementsByTagName("junction");
		for (int i = 0; i < nodesNodeList.getLength(); i++){
			Node node = nodesNodeList.item(i);
			if (node.getNodeType() == Node.ELEMENT_NODE){
				Element junctionElement = (Element) node;
				final String uuid = junctionElement.getAttribute("uuid");
				final String calleruuid = junctionElement.getAttribute("caller");
				final String chra = junctionElement.getAttribute("chra");
				final String chrb = junctionElement.getAttribute("chrb");
				final Long posa = Long.parseLong(junctionElement.getAttribute("posa"));
				final Long posb = Long.parseLong(junctionElement.getAttribute("posb"));
				Orientation oa = Orientation.valueOf(junctionElement.getAttribute("oapos"));
				Orientation ob = Orientation.valueOf(junctionElement.getAttribute("obpos"));
				Integer ebases = Integer.valueOf(junctionElement.getAttribute("ebases"));
				final int familyInt = junctionElement.hasAttribute("familyInt") ? Integer.parseInt(junctionElement.getAttribute("familyInt")) : 0;
				final int familyCount = junctionElement.hasAttribute("familyCount") ? Integer.parseInt(junctionElement.getAttribute("familyCount")) : 0;
				Caller thisCaller = SVGraph.callers.get(calleruuid);
				GraphJunction theJunction = new GraphJunction(chra, posa, chrb, posb, oa, ob, ebases, thisCaller, uuid);
				theJunction.setFamilyCount(familyCount);
				theJunction.setFamilyInt(familyInt);
				// append the partners
				NodeList partnersNodeList = junctionElement.getElementsByTagName("partner");
				for (int p = 0; p < partnersNodeList.getLength(); p++){
					Node pnode = partnersNodeList.item(p);
					if (pnode.getNodeType() == Node.ELEMENT_NODE){
						Element pElement = (Element) pnode;
						theJunction.addPartner(pElement.getAttribute("id"));
						theJunction.addPartnerCaller(pElement.getAttribute("caller"));
					}
				}

				// add the junction
				theGraph.addJunction(theJunction);
			}
		}

		return theGraph;

	}
	// =====================
	public void printCallers(){
		for (Caller c : callers.values()){
			System.out.println("Name: " + c.getName() + " Params: " + c.getParams());
		}
	}

	// =====================
	public boolean containsCaller(Caller c){
		return callers.values().contains(c);
	}

	// =====================
	public boolean containsCallerUUID(String s){
		return callers.containsKey(s);
	}

	// =====================
	public void sort(){
		// totally clear old sorting
		aSortedJunctions.clear();
		bSortedJunctions.clear();
		for (GraphJunction j : junctions.values()){
			if (!aSortedJunctions.containsKey(j.chra())){ aSortedJunctions.put(j.chra(), new ArrayList<String>()); }
			if (!bSortedJunctions.containsKey(j.chrb())){ bSortedJunctions.put(j.chrb(), new ArrayList<String>()); }
			aSortedJunctions.get(j.chra()).add(j.getUUID());
			bSortedJunctions.get(j.chrb()).add(j.getUUID());
		}
		for (List<String> cj : aSortedJunctions.values()){
			Collections.sort(cj, aComparator);
		}
		for (List<String> cj : bSortedJunctions.values()){
			Collections.sort(cj, bComparator);
		}
		sorted = true;
	}


	// =====================
	/**
	 * 		Annotates each GraphJunction as to weather or not it hits other junctions.  We have already separated the a and b sides of the junctions.
	 * 		
	 * 		sortedJunctions
	 * 			junctions are sorted by posl (the lower bound of the position
	 * 		j1	l===============h
	 * 		j2	 l====h
	 * 		j3	 l================h
	 * 		j4	               l=====h
	 * 		j5	                       l====h
	 * 		j6							l=====================h
	 * 
	 * 	  test	    l=====h
	 *    test	         l====h
	 *    
	 *    	We find all records such that test.l is less than q.l and stop looking when test.l is greater than q.h.
	 *     	Since we iterate over all records, this ensures that we will match all by all.
	 */
	public void matchUp(){
		System.err.println("Beginning matchup of " + junctions.size() + " junctions");
		// TODO; make this multithreading
		ExecutorService pool = Executors.newFixedThreadPool(Settings.threadCount);
		if (!sorted){ this.sort(); }
		for (List<String> gjl : aSortedJunctions.values()){
			pool.execute(new MergerWorker(gjl));
			//			System.err.println("Matching up " + gjl.size() + " records");
			//			for (String js : gjl){
			//				GraphJunction j = junctions.get(js);
			//				for (GraphJunction p : getMatchingJunctions(j)){
			//					j.addPartner(p.getUUID());
			//					j.addPartnerCaller(p.getCallerUUID());
			//					p.addPartner(j.getUUID());
			//					p.addPartnerCaller(j.getCallerUUID());
			//				}
			//			}
		}
		pool.shutdown();
		while (!pool.isTerminated()) { 
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {} 
		}
	}

	public void clearMatches(){
		for (GraphJunction j : junctions.values()){
			j.partners.clear();
		}
	}

	// =====================
	public void matchDown() throws Exception{
		if (!sorted){ this.sort(); }
		Map<String, List<GraphJunction>> lJunctionsByCaller = new HashMap<String, List<GraphJunction>>();
		for (GraphJunction j : getJunctions()){
			final String caller = j.getCallerUUID();
			if (!lJunctionsByCaller.containsKey(caller)){ lJunctionsByCaller.put(caller, new ArrayList<GraphJunction>()); }
			lJunctionsByCaller.get(caller).add(j);
		}
		ExecutorService pool = Executors.newFixedThreadPool(Settings.threadCount);
		for (String ck : lJunctionsByCaller.keySet()){
			pool.execute(new MatchDownWorker(lJunctionsByCaller.get(ck), ck));

			//			System.out.println("Matching down junctions from caller; " + SVGraph.callers.get(ck).getName());
			//			List<GraphJunction> callerJunctionsASorted = new ArrayList<GraphJunction>();
			//			List<GraphJunction> callerJunctionsBSorted = new ArrayList<GraphJunction>();
			//			List<GraphJunction> junctions = lJunctionsByCaller.get(ck);
			//			callerJunctionsASorted.addAll(lJunctionsByCaller.get(ck));
			//			callerJunctionsBSorted.addAll(lJunctionsByCaller.get(ck));
			//			Collections.sort(callerJunctionsASorted, GraphJunction.aComparator);
			//			Collections.sort(callerJunctionsBSorted, GraphJunction.bComparator);
			//			Set<GraphJunction> seenSet = new HashSet<GraphJunction>();
			//			List<GraphJunction> addingJunctions = new ArrayList<GraphJunction>();
			//			List<GraphJunction> toRemove = new ArrayList<GraphJunction>();
			//			for (GraphJunction j : junctions){
			//				// find the family for this junction
			//				if (seenSet.contains(j)){ continue; }
			//				Set<GraphJunction> newSet = new HashSet<GraphJunction>();
			//				newSet.add(j);
			//				_popMatchDownGJSet(newSet, j, callerJunctionsASorted, callerJunctionsBSorted);
			//				seenSet.addAll(newSet);
			//				if (newSet.size() > 1){
			//					toRemove.addAll(newSet); // we need to remove all of this family of junctions
			//					System.out.println("Matchdown removing " + newSet.size() + " junctions.");
			//					GraphJunction newJunction = _combineGraphJunctions(newSet);
			//					addingJunctions.add(newJunction);
			//				}
			//			}
			//			for (GraphJunction rj : toRemove){
			//				SVGraph.junctions.remove(rj.getUUID());
			//			}
			//			for (GraphJunction aj : addingJunctions){
			//				this.addJunction(aj);
			//			}
		}
		pool.shutdown();
		while (!pool.isTerminated()) { 
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {} 
		}
	}

	private static GraphJunction _combineGraphJunctions(Set<GraphJunction> newSet) {
		List<GraphJunction> oldJunctions = new ArrayList<GraphJunction>();
		oldJunctions.addAll(newSet);
		// find the averages for the set and generate a new graph junction
		GraphJunction firstJunction = oldJunctions.get(0);
		final String chra = firstJunction.chra();
		final String chrb = firstJunction.chrb();
		final Caller caller = SVGraph.callers.get(firstJunction.getCallerUUID());
		final Orientation oa = firstJunction.oapos();
		final Orientation ob = firstJunction.obpos();
		final boolean samechr = chra.equals(chrb); // we only need to calculate this once
		long posa = 0;
		long posb = 0;
		int ebases = 0;
		for (GraphJunction oj : oldJunctions){
			if (samechr){
				if (oj.posa() < oj.posb()){
					posa = posa + oj.posa();
					posb = posb + oj.posb();
				} else {
					posa = posa + oj.posb();
					posb = posa + oj.posa();
				}
			} else {
				if (chra.equals(oj.chra()) && chrb.equals(oj.chrb())){
					posa = posa + oj.posa();
					posb = posb + oj.posb();
				} else if ( chrb.equals(oj.chra()) && chra.equals(oj.chrb()) ){
					posa = posa + oj.posb();
					posb = posa + oj.posa();
				} else {
					// sanity check this should never happen
					System.err.println("Something really bad just happened, I tried to merge a junciton from a different cluster");
					System.err.println(" arch chra: " + chra + " chrb: " + chrb + " new chra: " + oj.chra() + " chrb: " + oj.chrb());
					continue;
				}
			}
			ebases = ebases + oj.ebases();
		}
		return new GraphJunction(chra, posa / oldJunctions.size(), chrb, posb / oldJunctions.size(), oa, ob, ebases / oldJunctions.size(), caller );
	}

	private static void _popMatchDownGJSet(Set<GraphJunction> thisSet, GraphJunction j, List<GraphJunction> callerJunctionsASorted, List<GraphJunction> callerJunctionsBSorted){
		int aaHint = getAHintGJ(callerJunctionsASorted, j.chra(), j.posa(), 0, callerJunctionsASorted.size() - 1);
		for (int i = aaHint; i < callerJunctionsASorted.size(); i++){
			GraphJunction tj = callerJunctionsASorted.get(i);
			if ((tj.posa() - j.posa() > 10) || (!tj.chra().equals(j.chra()))){ break; }
			if (thisSet.contains(tj)){ continue; }
			if (localIntersects(j, tj)){
				thisSet.add(tj);
				_popMatchDownGJSet(thisSet, tj, callerJunctionsASorted, callerJunctionsBSorted);
			}
		}
		int abHint = getBHintGJ(callerJunctionsBSorted, j.chra(), j.posa(), 0, callerJunctionsBSorted.size() - 1);
		for (int i = abHint; i < callerJunctionsBSorted.size(); i++){
			GraphJunction tj = callerJunctionsBSorted.get(i);
			if ((tj.posb() - j.posa() > 10) || (!tj.chrb().equals(j.chra()))){ break; }
			if (thisSet.contains(tj)){ continue; }
			if (localIntersects(j, tj)){
				thisSet.add(tj);
				_popMatchDownGJSet(thisSet, tj, callerJunctionsASorted, callerJunctionsBSorted);
			}
		}
		int baHint = getAHintGJ(callerJunctionsASorted, j.chrb(), j.posb(), 0, callerJunctionsASorted.size() - 1);
		for (int i = baHint; i < callerJunctionsASorted.size(); i++){
			GraphJunction tj = callerJunctionsASorted.get(i);
			if ((tj.posa() - j.posb() > 10) || (!tj.chra().equals(j.chrb()))){ break; }
			if (thisSet.contains(tj)){ continue; }
			if (localIntersects(j, tj)){
				thisSet.add(tj);
				_popMatchDownGJSet(thisSet, tj, callerJunctionsASorted, callerJunctionsBSorted);
			}
		}
		int bbHint = getBHintGJ(callerJunctionsBSorted, j.chrb(), j.posb(), 0, callerJunctionsBSorted.size() - 1);
		for (int i = bbHint; i < callerJunctionsBSorted.size(); i++){
			GraphJunction tj = callerJunctionsBSorted.get(i);
			if ((tj.posb() - j.posb() > 10) || (!tj.chrb().equals(j.chrb()))){ break; }
			if (thisSet.contains(tj)){ continue; }
			if (localIntersects(j, tj)){
				thisSet.add(tj);
				_popMatchDownGJSet(thisSet, tj, callerJunctionsASorted, callerJunctionsBSorted);
			}
		}
	}

	private static boolean localIntersects(GraphJunction a, GraphJunction b){
		if (a.chra().equals(b.chra()) && a.chrb().equals(b.chrb())){
			if ((Math.abs(a.posa() - b.posa()) < 10) && (Math.abs(a.posb() - b.posb()) < 10)){
				return true;
			} else {
				return false;
			}
		} else if (a.chra().equals(b.chrb()) && a.chrb().equals(b.chra())){
			if ((Math.abs(a.posa() - b.posb()) < 10) && (Math.abs(a.posb() - b.posa()) < 10)){
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	// =====================
	private static int getAHintGJ(List<GraphJunction> gl, String chr, long pos, int min, int max){
		if (max < min){ return max > 0 ? max : 0; }
		int mid = (max + min) / 2;
		final GraphJunction val = gl.get(mid);
		final int chrComp = chr.compareTo(val.chra());
		if (chrComp < 0){
			return getAHintGJ(gl, chr, pos, min, mid -1);
		} else if (chrComp > 0){
			return getAHintGJ(gl, chr, pos, mid + 1, max);
		} else {
			// compare the position
			if (pos > val.posa()){ return getAHintGJ(gl, chr, pos, mid + 1, max); }
			else { return getAHintGJ(gl, chr, pos, min, mid -1); }
		}
	}

	// =====================
	private static int getBHintGJ(List<GraphJunction> gl, String chr, long pos, int min, int max){
		if (max < min){ return max > 0 ? max : 0; }
		int mid = (max + min) / 2;
		final GraphJunction val = gl.get(mid);
		final int chrComp = chr.compareTo(val.chra());
		if (chrComp < 0){
			return getBHintGJ(gl, chr, pos, min, mid -1);
		} else if (chrComp > 0){
			return getBHintGJ(gl, chr, pos, mid + 1, max);
		} else {
			// compare the position
			if (pos > val.posa()){ return getBHintGJ(gl, chr, pos, mid + 1, max); }
			else { return getBHintGJ(gl, chr, pos, min, mid -1); }
		}
	}

	// =====================
	/**
	 * Return the position of the 
	 * @param chr
	 * @param pos
	 * @return
	 */
	private static int getAHint(List<String> gl, long pos, int min, int max){
		if (max < min){ return max > 0 ? max : 0; }
		int mid = (max + min) / 2;
		final String vals = gl.get(mid);
		final GraphJunction val = junctions.get(vals);
		long tpos = val.getPosal();
		if (tpos > pos){ return getAHint(gl, pos, min, mid -1); }
		else if (tpos < pos){ return getAHint(gl, pos, mid + 1, max); }
		else {
			// they are equal, we back track until pos is less than test
			int i = mid - 1;
			while(i > 0){
				final GraphJunction bt = junctions.get(gl.get(i));
				if (bt.getPosal() < pos){ return i; }
				else { i--; }
			} ;
			return i > 0 ? i : 0;
		}
	}

	// =====================
	private static int getBHint(List<String> gl, long pos, int min, int max){
		if (max < min){ return max > 0 ? max : 0; }
		int mid = (max + min) / 2;
		final String vals = gl.get(mid);
		final GraphJunction val = junctions.get(vals);
		long tpos = val.getPosbl();
		if (tpos > pos){ return getBHint(gl, pos, min, mid - 1); }
		else if (tpos < pos){ return getBHint(gl, pos, mid + 1, max); }
		else {
			// they are equal, we back track until pos is less than test
			int i = mid - 1;
			while (i > 0){
				final GraphJunction bt = junctions.get(gl.get(i));
				if (bt.getPosbl() < pos){ return i; }
				else { i--; }
			} ;
			return i > 0 ? i : 0;
		}
	}

	// =====================
	public Set<GraphJunction> getMatchingJunctions(GraphJunction j){
		if (!sorted){ this.sort(); }
		Map<String, List<GraphJunction>> callerMatches = new HashMap<String, List<GraphJunction>>();
		Set<GraphJunction> resultSet = new HashSet<GraphJunction>();
		final String achr = j.chra();
		final String bchr = j.chrb(); 
		final List<String> aachrList = aSortedJunctions.get(achr);
		final List<String> abchrList = bSortedJunctions.get(achr);
		final List<String> bachrList = aSortedJunctions.get(bchr);
		final List<String> bbchrList = bSortedJunctions.get(bchr);
		if (aachrList != null){
			int aaHint = SVGraph.getAHint(aachrList, j.getPosal(), 0, aachrList.size() - 1);
			// process the aa junctions
			for (int aa = aaHint; aa < aachrList.size(); aa++){
				final String jts = aachrList.get(aa);
				final GraphJunction jt = junctions.get(jts);
				if (jt.getCallerUUID().equals(j.getCallerUUID())){ continue; } // don't add self or near self.
				if (!callerMatches.containsKey(jt.getCallerUUID())){ callerMatches.put(jt.getCallerUUID(), new ArrayList<GraphJunction>()); }
				if (j.intersects(jt)){
					callerMatches.get(jt.getCallerUUID()).add(jt);
				} else if (j.getPosal() > jt.getPosah()){
					break;
				}
			}
		}

		if (abchrList != null){
			int abHint = SVGraph.getBHint(abchrList, j.getPosal(), 0, abchrList.size() - 1);
			// process the ab junctions
			for (int ab = abHint; ab < abchrList.size(); ab++){
				final String jts = abchrList.get(ab);
				final GraphJunction jt = junctions.get(jts);
				if (jt.getCallerUUID().equals(j.getCallerUUID())){ continue; } // don't add self or near self.
				if (!callerMatches.containsKey(jt.getCallerUUID())){ callerMatches.put(jt.getCallerUUID(), new ArrayList<GraphJunction>()); }
				if (j.intersects(jt)){
					callerMatches.get(jt.getCallerUUID()).add(jt);
				} else if (j.getPosal() > jt.getPosbh()){
					break;
				}
			}
		}

		if (bachrList != null){
			int baHint = SVGraph.getAHint(bachrList, j.getPosbl(), 0, bachrList.size() - 1);
			// process the ba junctions
			for (int ba = baHint; ba < bachrList.size(); ba++){
				final String jts = bachrList.get(ba);
				final GraphJunction jt = junctions.get(jts);
				if (jt.getCallerUUID().equals(j.getCallerUUID())){ continue; } // don't add self or near self.
				if (!callerMatches.containsKey(jt.getCallerUUID())){ callerMatches.put(jt.getCallerUUID(), new ArrayList<GraphJunction>()); }
				if (j.intersects(jt)){
					callerMatches.get(jt.getCallerUUID()).add(jt);
				} else if (j.getPosbl() > jt.getPosah()){
					break;
				}
			}
		}

		if (bbchrList != null){
			int bbHint = SVGraph.getBHint(bbchrList, j.getPosbl(), 0, bbchrList.size() - 1);
			// process the bb junctions
			for (int bb = bbHint; bb < bbchrList.size(); bb++){
				final String jts = bbchrList.get(bb);
				final GraphJunction jt = junctions.get(jts);
				if (jt.getCallerUUID().equals(j.getCallerUUID())){ continue; } // don't add self or near self.
				if (!callerMatches.containsKey(jt.getCallerUUID())){ callerMatches.put(jt.getCallerUUID(), new ArrayList<GraphJunction>()); }
				if (j.intersects(jt)){
					callerMatches.get(jt.getCallerUUID()).add(jt);
				} else if (j.getPosbl() > jt.getPosbh()){
					break;
				}
			}
		}

		for (String ck : callerMatches.keySet()){
			// we only allow at most two matches per call set, and only one if they are full hits
			List<GraphJunction> gjList = callerMatches.get(ck);
			GraphJunction fullJunction = null;
			GraphJunction leftJunction = null;
			GraphJunction rightJunction = null;

			for (GraphJunction cj : gjList){
				HitType h = j.hitType(cj);
				if (h.equals(HitType.Full)){ 
					fullJunction = cj;
					break; // don't need to go farther, we are breaking
				} else if (h.equals(HitType.Left)){
					leftJunction = cj;
				} else if (h.equals(HitType.Right)){
					rightJunction = cj;
				}
			}
			if (fullJunction != null){
				resultSet.add(fullJunction);
			} else {
				if (leftJunction != null){ resultSet.add(leftJunction); }
				if (rightJunction != null){ resultSet.add(rightJunction); }
			}

		}

		return resultSet;
	}

	// =====================
	public int nodeCount() {
		return junctions.size();
	}

	// =====================
	public int edgeCount() {
		int count = 0;
		for (GraphJunction g : junctions.values()){
			count += g.partnerCount();
		}
		return count / 2;
	}

	// =====================
	public Set<String> getPartnerCallers(String key){
		final GraphJunction j = junctions.get(key);
		return getPartnerCallers(j);
	}

	// =====================
	public Set<String> getPartnerCallers(GraphJunction j){
		Set<String> result = new HashSet<String>();
		if (j == null) return result;
		for (String p : j.partners){
			final GraphJunction t = junctions.get(p);
			result.add(t.getCallerUUID());
		}
		return result;
	}

	// =====================
	public Caller getCaller(String key) {
		return callers.get(key);
	}

	// =====================
	public Caller getCallerByName(String key) {
		for (Caller c : callers.values()){
			if (c.getName().equals(key)){ return c; }
		}
		return null;
	}

	// =====================
	public Map<String, List<GraphJunction>> getCallerData() {
		Map<String, List<GraphJunction>> result = new HashMap<String, List<GraphJunction>>();
		for (GraphJunction j : junctions.values()){
			if (!result.containsKey(j.getCallerUUID())){ result.put(j.getCallerUUID(), new ArrayList<GraphJunction>()); }
			result.get(j.getCallerUUID()).add(j);
		}
		return result;
	}

	// =====================
	public List<GraphJunction> getJunctions(){
		List<GraphJunction> result = new ArrayList<GraphJunction>(junctions.size());
		result.addAll(junctions.values());
		return result;
	}

}
