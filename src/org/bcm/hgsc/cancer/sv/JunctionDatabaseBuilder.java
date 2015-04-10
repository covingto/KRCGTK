package org.bcm.hgsc.cancer.sv;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.List;
import java.util.Map;

import org.json.simple.JSONValue;

public class JunctionDatabaseBuilder {
	
	private static String sCreateEnds = "CREATE TABLE IF NOT EXISTS ENDS (k INTEGER PRIMARY KEY ASC, chr TEXT, lpos INT, pos INT, hpos INT, orientation TEXT)";
	private static String sCreateJunctions = "CREATE TABLE IF NOT EXISTS JUNCTIONS (k INTEGER PRIMARY KEY ASC, lend INT, rend INT, caller INT, ebases INT)";
	private static String sCreateCallers = "CREATE TABLE IF NOT EXISTS CALLERS (k INTEGER PRIMARY KEY ASC, name TEXT, params TEXT)";
	private static String psAddCaller = "INSERT INTO CALLERS (name, params) VALUES (?, ?)";
	private static String psAddEnd = "INSERT INTO ENDS (chr, lpos, pos, hpos, orientation) VALUES (?, ?, ?, ?, ?)";
	private static String psAddJunction = "INSERT INTO JUNCTIONS (lend, rend, caller, ebases) VALUES (?, ?, ?, ?)";
	private static String indexEndsAll = "CREATE INDEX ENDINDEX ON ENDS (chr, pos)";
	private static String indexEndsChr = "CREATE INDEX CHRINDEX ON ENDS (chr)";
	private static String indexEndsPos = "CREATE INDEX POSINDEX ON ENDS (pos)";
	private static String indexEndsH = "CREATE INDEX POSINDEXH ON ENDS (hpos)";
	private static String indexEndsL = "CREATE INDEX POSINDEXL ON ENDS (lpos)";
	private static String indexJunctionLkeys = "CREATE INDEX JUNCTIONLKEYS ON JUNCTIONS (lend)";
	private static String indexJunctionRkeys = "CREATE INDEX JUNCTIONRKEYS ON JUNCTIONS (rend)";
	private static String indexJunctionsCaller = "CREATE INDEX JUNCTIONCALLERS ON JUNCTIONS (caller)";
//	private static String createLeftView = "CREATE VIEW LEFTJUNCTIONS AS SELECT ENDS.k AS KEYL, ENDS.chr AS CHRL, ENDS.pos AS POSL, ENDS.orientation AS ORL FROM JUNCTIONS, ENDS ON ENDS.k = JUNCTIONS.lend";
//	private static String createRightView = "CREATE VIEW RIGHTJUNCTIONS AS SELECT ENDS.k AS KEYR, ENDS.chr AS CHRR, ENDS.pos AS POSR, ENDS.orientation AS ORR FROM JUNCTIONS, ENDS ON ENDS.k = JUNCTIONS.rend";
	private static String createJunctionTable = "CREATE TABLE FULLJUNCTIONS AS SELECT DISTINCT " +
			"CHRL, LPOSL, POSL, HPOSL, CHRR, LPOSR, POSR, HPOSR, CALLERS.name AS CALLER, CALLERS.k as CALLERKEY, JUNCTIONS.k AS JKEY, KEYL, KEYR " +
			"FROM JUNCTIONS, " +
			"(SELECT ENDS.k AS KEYL, ENDS.chr AS CHRL, ENDS.lpos AS LPOSL, ENDS.hpos AS HPOSL, ENDS.pos AS POSL, ENDS.orientation AS ORL FROM JUNCTIONS, ENDS ON ENDS.k = JUNCTIONS.lend) AS LJUNCTIONS, " +
			"(SELECT ENDS.k AS KEYR, ENDS.chr AS CHRR, ENDS.lpos AS LPOSR, ENDS.hpos AS HPOSR, ENDS.pos AS POSR, ENDS.orientation AS ORR FROM JUNCTIONS, ENDS ON ENDS.k = JUNCTIONS.rend) AS RJUNCTIONS, " +
			"CALLERS " +
			"ON CALLERS.k = JUNCTIONS.caller AND JUNCTIONS.lend = LJUNCTIONS.KEYL AND JUNCTIONS.rend = RJUNCTIONS.KEYR";
	private static String indexFullJunctionLkey = "CREATE INDEX FULLKEYL ON FULLJUNCTIONS (KEYL)";
	private static String indexFullJunctionRkey = "CREATE INDEX FULLKEYR ON FULLJUNCTIONS (KEYR)";
	private static String indexFullJunctionLloc = "CREATE INDEX FULLLLOC ON FULLJUNCTIONS (CHRL, LPOSL, POSL, HPOSL)";
	private static String indexFullJunctionRloc = "CREATE INDEX FULLRLOC ON FULLJUNCTIONS (CHRR, LPOSR, POSR, HPOSR)";
	private static String indexFullJunctionLc = "CREATE INDEX FULLLC ON FULLJUNCTIONS (CHRL)";
	private static String indexFullJunctionRc = "CREATE INDEX FULLRC ON FULLJUNCTIONS (CHRR)";
	private static String indexFullJunctionLp = "CREATE INDEX FULLLP ON FULLJUNCTIONS (POSL)";
	private static String indexFullJunctionRp = "CREATE INDEX FULLRP ON FULLJUNCTIONS (POSR)";
	private static String indexFullJunctionLkloc = "CREATE INDEX FULLLKLOC ON FULLJUNCTIONS (KEYL, CHRL, POSL)";
	private static String indexFullJunctionRkloc = "CREATE INDEX FULLRKLOC ON FULLJUNCTIONS (KEYR, CHRR, POSR)";
	private static String createCrossHits = "CREATE TABLE CROSSHITS AS " +
			"SELECT P.CALLERKEY AS PCALLER, H.CALLERKEY AS HCALLER, P.JKEY AS PKEY, H.JKEY AS HKEY, " +
			"(CASE " +
				"WHEN " +
				"((P.CHRL == H.CHRL AND P.LPOSL < H.HPOSL AND P.HPOSL > H.LPOSL) OR " +
				"(P.CHRL == H.CHRR AND P.LPOSL < H.HPOSR AND P.HPOSL > H.LPOSR)) " +
				"THEN 1 " +
				"ELSE 0 END) AS LEFTHIT, " +
			"(CASE " +
				"WHEN " +
				"((P.CHRR == H.CHRL AND P.LPOSR < H.HPOSL AND P.HPOSR > H.LPOSL) OR " +
				"(P.CHRR == H.CHRR AND P.LPOSR < H.HPOSR AND P.HPOSR > H.LPOSR)) " +
				"THEN 1 " +
				"ELSE 0 END) AS RIGHTHIT " +
			"FROM FULLJUNCTIONS AS P LEFT JOIN FULLJUNCTIONS AS H " +
			"ON (P.CALLERKEY < H.CALLERKEY) AND " +
			"(" +
			"(P.CHRL == H.CHRL AND P.LPOSL < H.HPOSL AND P.HPOSL > H.LPOSL) OR " +
			"(P.CHRL == H.CHRR AND P.LPOSL < H.HPOSR AND P.HPOSL > H.LPOSR) OR " +
			"(P.CHRR == H.CHRL AND P.LPOSR < H.HPOSL AND P.HPOSR > H.LPOSL) OR " +
			"(P.CHRR == H.CHRR AND P.LPOSR < H.HPOSR AND P.HPOSR > H.LPOSR)" +
			")";
	private static String indexCrossPkey = "CREATE INDEX CROSSPKEY ON CROSSHITS (PKEY)";
	private static String indexCrossHkey = "CREATE INDEX CROSSHKEY ON CROSSHITS (HKEY)";
	private static String indexCrossHtypeL = "CREATE INDEX CROSSHTYPEL ON CROSSHITS (LEFTHIT)";
	private static String indexCrossHtypeR = "CREATE INDEX CROSSHTYPER ON CROSSHITS (RIGHTHIT)";

	public static void main(String[] args) throws Exception {
		//Options options = new Options();
		//Parser parser = new BasicParser();
		
		//options.addOption("c", false, "Configuration file indicates database creation, move to creation stream");
		
		Connection c = null;
	    
		if (args.length < 1){
			System.err.println("Useage: JunctionDatabaseBuilder.jar conffile.conf");
			System.exit(11);
		}
	    File conffile = new File(args[0]);
	    
	    try {
	      Class.forName("org.sqlite.JDBC");
	      System.out.println("Opened database successfully");
	      FileReader fr = new FileReader(conffile);
	      @SuppressWarnings("unchecked")
	      Map<String, Object> conf = (Map<String, Object>) JSONValue.parse(fr);
	      fr.close();
	      String dbname = (String) conf.get("dbname");
	      c = DriverManager.getConnection("jdbc:sqlite:" + dbname);
	      c.setAutoCommit(false); // theoretical performance boost
	      
	      @SuppressWarnings("unchecked")
	      List<Object> settings = (List<Object>) conf.get("settings");
	      for (Object s : settings){
	    	  @SuppressWarnings("unchecked")
			Map<String, Object> m = (Map<String, Object>) s;
	    	  appendData(c, (String) m.get("source"), (String) m.get("caller"), (String) m.get("conf"), (Long) m.get("buffer"), (String) m.get("platform"));
	      }
	      
	      c.prepareStatement(indexEndsAll).executeUpdate();
	      c.prepareStatement(indexEndsChr).executeUpdate();
	      c.prepareStatement(indexEndsPos).executeUpdate();
	      c.prepareStatement(indexEndsH).executeUpdate();
	      c.prepareStatement(indexEndsL).executeUpdate();
	      c.prepareStatement(indexJunctionLkeys).executeUpdate();
	      c.prepareStatement(indexJunctionRkeys).executeUpdate();
	      c.prepareStatement(indexJunctionsCaller).executeUpdate();
	      
	      System.out.println("Creating junction table");
	      c.prepareStatement(createJunctionTable).executeUpdate();
	      
	      c.prepareStatement(indexFullJunctionLkey).executeUpdate();
	      c.prepareStatement(indexFullJunctionRkey).executeUpdate();
	      c.prepareStatement(indexFullJunctionLloc).executeUpdate();
	      c.prepareStatement(indexFullJunctionRloc).executeUpdate();
	      c.prepareStatement(indexFullJunctionLc).executeUpdate();
	      c.prepareStatement(indexFullJunctionRc).executeUpdate();
	      c.prepareStatement(indexFullJunctionLp).executeUpdate();
	      c.prepareStatement(indexFullJunctionRp).executeUpdate();
	      c.prepareStatement(indexFullJunctionLkloc).executeUpdate();
	      c.prepareStatement(indexFullJunctionRkloc).executeUpdate();
	      
	      System.out.println("Creating crosshits");
	      c.prepareStatement(createCrossHits).executeUpdate();
	      
	      c.prepareStatement(indexCrossPkey).executeUpdate();
	      c.prepareStatement(indexCrossHkey).executeUpdate();
	      c.prepareStatement(indexCrossHtypeL).executeUpdate();
	      c.prepareStatement(indexCrossHtypeR).executeUpdate();
	      
	      System.out.println("Committing data to database");
	      
	      c.commit();
	      c.close();
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	    System.out.println("Database generated successfully.");
	}

	private static void appendData(Connection c, String source, String caller, String conf, long buffer, String platform) throws Exception {
		// create tables if they do not exist
		System.out.println("Parsing data;");
		System.out.println("Source: " + source);
		System.out.println("Caller: " + caller);
		System.out.println("Conf: " + conf);
		try {
			c.prepareStatement(sCreateCallers).executeUpdate();
		} catch (Exception eCaller) {
			System.err.println("Exception in creating callers table");
			System.err.println(sCreateCallers);
			eCaller.printStackTrace();
		}
		try {
			c.prepareStatement(sCreateEnds).executeUpdate();
		} catch (Exception eEnds) {
			System.err.println("Exception in creating ends table");
			System.err.println(sCreateEnds);
			eEnds.printStackTrace();
		}
		try {
			c.prepareStatement(sCreateJunctions).executeUpdate();
		} catch (Exception eJunctions) {
			System.err.println("Exception in creating junctions table");
			System.err.println(sCreateJunctions);
			eJunctions.printStackTrace();
		}
		
		// add the caller info to the database
		try {
			PreparedStatement cs = c.prepareStatement(psAddCaller, Statement.RETURN_GENERATED_KEYS);
			PreparedStatement ae = c.prepareStatement(psAddEnd, Statement.RETURN_GENERATED_KEYS);
			PreparedStatement aj = c.prepareStatement(psAddJunction, Statement.RETURN_GENERATED_KEYS);
			
			cs.setString(1, caller);
			cs.setString(2, conf);
			cs.executeUpdate();
			ResultSet callerIDSet = cs.getGeneratedKeys();
			callerIDSet.next();
			int callerInt = callerIDSet.getInt(1);
			
			// parse through the input file and add the records
			Caller theCaller = new Caller(caller, conf, 1000, platform);
			try {
				JunctionSourceReader jsr;
				if (caller.equals("BreakDown")){
					jsr = new BreakDownJunctionReader(source, theCaller);
				} else if (caller.equals("BreakDancer")){
					jsr = new BreakDancerReader(source, theCaller);
				} else if (caller.equals("PInDel")){
					jsr = new PInDelJunctionReader(source, theCaller);
				} else if (caller.equals("GATK")) {
					jsr = new VCFJunctionReader(source, theCaller);
				} else if (caller.equals("HoneyTails")) {
					jsr = new HoneyTailsJunctionReader(source, theCaller);
				} else if (caller.equals("HoneySpots")) {
					jsr = new HoneySpotsJunctionReader(source, theCaller);
				} else {
					jsr = new GenericJunctionReader(source, theCaller);
				}
				int nrecs = 0;
				while (jsr.hasNext()){
					final JRecord jr = jsr.next();
					
					ae.setString(1, jr.chra());
					ae.setFloat(2, jr.posa() - buffer);
					ae.setFloat(3, jr.posa());
					ae.setFloat(4, jr.posa() + buffer);
					ae.setString(5, jr.oapos().toString());
					ae.executeUpdate();
					final int lkey = key(ae.getGeneratedKeys());
					
					ae.setString(1, jr.chrb());
					ae.setFloat(2, jr.posb() - buffer);
					ae.setFloat(3, jr.posb());
					ae.setFloat(4, jr.posb() + buffer);
					ae.setString(5, jr.obpos().toString());
					ae.executeUpdate();
					final int rkey = key(ae.getGeneratedKeys());
					
					aj.setInt(1, lkey);
					aj.setInt(2, rkey);
					aj.setInt(3, callerInt);
					aj.setInt(4, jr.ebases());
					aj.executeUpdate();
					nrecs++;
				}
				jsr.close();
				System.out.println("Wrote " + nrecs + " to database");
			} catch (IOException readExcept){
				System.err.println("Error reading source file: " + source);
				readExcept.printStackTrace();
				throw readExcept;
			} catch (SQLException sqlExcept) {
				System.err.println("SQL Exception on table add");
				sqlExcept.printStackTrace();
				throw sqlExcept;
			}
		} catch (SQLException eCallerAdd){
			System.err.println("Failure to add caller information");
			System.err.println(psAddCaller + ", " + caller + ", " + conf);
			eCallerAdd.printStackTrace();
			throw eCallerAdd;
		}
		
	}
	
	private static int key(ResultSet k) throws SQLException{
		k.next();
		return k.getInt(1);
	}
}
