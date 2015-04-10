package org.bcm.hgsc.cancer;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.cancer.utils.Mafstruct;

public class GermlineFilter {
	
	private static final List<String> indels = Arrays.asList(new String[] {"INS", "DEL", "COMPLEX"});
	private enum compRes { LOWER, HIGHER, SNVHIT, INDELHIT }
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// args should be a germline and a somatic maf file
		
		// this runs more like a script than anything we just use the power of Java's collections
		try {
			 /*
			CommandLine commandLine;
			@SuppressWarnings("static-access")
			Option optSomaticPass = OptionBuilder.withArgName("somaticPassOut")
					.hasArg()
					.withDescription("Output file name for the passing variants.")
					.create("somatic-pass");
			@SuppressWarnings("static-access")
			Option optSomaticFail = OptionBuilder.withArgName("somaticFailOut")
					.hasArg()
					.withDescription("Output file name for the failing variants.")
					.create("somatic-fail");
			
			Options options = new Options();
			CommandLineParser parser = new BasicParser();
			
			options.addOption(optSomaticFail);
			options.addOption(optSomaticPass);
			
			commandLine = parser.parse(options, args);
			*/
			
			System.out.println("Running germline removal with the following args;");
			System.out.println(StringUtils.join(args, " "));
			
			String germlinefile = args[0];
			String somaticfile = args[1];
			String somaticPassOut = args[2];
			String somaticFlag = args[3];
			
			System.out.println("Generating maf structures.");
			Mafstruct germline = new Mafstruct(new File(germlinefile));
			Mafstruct somatic = new Mafstruct(new File(somaticfile));
			
			System.out.println("Flagging variants.");
			linearVarHit(germline, somatic);
			
			System.out.println("Writing to files.");
			somatic.writePassing(somaticPassOut);
			somatic.writeFailing(somaticFlag);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void germlineFlag(Mafstruct germline, Mafstruct somatic){
		Map<String, List<MAFData>> subjectHashG = new HashMap<String, List<MAFData>>();
		Map<String, List<MAFData>> subjectHashS = new HashMap<String, List<MAFData>>();
		
		// build the hashes
		System.out.println("Generating subject hashes");
		for (String s : germline.subjects()){
			subjectHashG.put(s, germline.subsetSubject(s));
		}
		for (String s : somatic.subjects()){
			subjectHashS.put(s, somatic.subsetSubject(s));
		}
		
		// flag variants germline
		System.out.println("Flagging variants.");
		int pass = 0;
		int selfgermline = 0;
		int othergermline = 0;
		for (MAFData v : somatic.getVariants()){
			// test a match against self
			if ( binaryVarHit( subjectHashG.get(v.get("Tumor_Sample_Barcode")), v )  < 0 ) {
				if ( binaryVarHit( germline, v) < 0){
					v.set("Failed_Reason", "Pass");
					pass += 1;
				} else {
					v.set("Failed_Reason", "OtherGermline");
					othergermline += 1;
				}
			} else {
				v.set("Failed_Reason", "SelfGermline");
				selfgermline += 1;
			}
		}
		System.out.println("Found " + pass + " passing variants");
		System.out.println("Found " + selfgermline + " self germline variants");
		System.out.println("Found " + othergermline + " other germline variants");
	}
	
	public static void linearVarHit(Mafstruct germline, Mafstruct somatic) {
		System.out.println("Sorting germline and somatic variants");
		germline.sort();
		somatic.sort();
		List<MAFData> gvars = germline.getVariants();
		List<MAFData> svars = somatic.getVariants();
		int i = 0; // the somatic counter
		int j = 0; // the germline counter
		for (i = 0; i < svars.size(); i++){
			j = _linearSearchGermline(svars.get(i), gvars, j);
		}
	}
	
	private static int _linearSearchGermline(MAFData var, List<MAFData> s, int j){
		int k = j;
		
		String vsub = var.get("Tumor_Sample_Barcode");
		
		int mink = k;
		boolean match = false;
		
		while(!match & k < s.size()) {
			MAFData sub = s.get(k);
			int intersect = intersects(sub, var);
			if (intersect < 0) {
				// the subject is before the variant in a list
				k++;
				mink = k;
			} else if (intersect > 0) {
				// the subject is after the variant in a list
				// this means we exceed the counter
				return mink;
			} else {
				// intersect hits
				// see if this is the same subject, if not 
				String ssub = sub.get("Tumor_Sample_Barcode");
				if (ssub.equals(vsub)){
					var.set("Failed_Reason", "SelfGermline");
					return mink;
				} else {
					var.set("Failed_Reason", "OtherGermline");
					// but we continue in case we find a self germline
					k++;
				}
			}
		} 
		
		return mink;
	}
	
	public static int binaryVarHit(List<MAFData> subject, MAFData variant) {
		return binaryVarHit(subject, variant, 0, subject.size());
	}
	
	public static int binaryVarHit(Mafstruct subject, MAFData variant) {
		return binaryVarHit(subject.getVariants(), variant, 0, subject.length());
	}
	
	public static int binaryVarHit(List<MAFData> subject, MAFData variant, int imin, int imax){
		
		if ( imax < imin ){ return -imin; }
		int imid = (imin + imax) / 2;
		MAFData s = subject.get(imid);
		
		// test if the key intersects the s
		int hit = intersects(s, variant);
		if ( hit > 0 ){
			return binaryVarHit(subject, variant, imin, imid - 1);
		} else if ( hit < 0 ){
			return binaryVarHit(subject, variant, imid + 1, imax);
		} else {
			// key has been found
			return imid;
		}
	}
	
	public static int intersects(MAFData s, MAFData v) {
		if ( indels.contains(s.get("Variant_Type")) ){
			return indelIntersect(s, v);
		} else {
			return snvIntersect(s, v);
		}
	}
	
	public static boolean between(int a, int b, int q){
		if (q <= b & q >= a) { return true; }
		return false;
	}
	
	public static int indelIntersect(MAFData s, MAFData v) {
		// the subject is an indel so we expand the footprint by INDELBUFFER
		int chrcomp = s.get("Chromosome").compareTo(v.get("Chromosome"));
		if (chrcomp == 0){
			Integer sstart = Integer.parseInt(s.get("Start_position")) - 5;
			Integer send = Integer.parseInt(s.get("End_position")) + 5;
			Integer vstart = Integer.parseInt(v.get("Start_position"));
			Integer vend = Integer.parseInt(v.get("End_position"));
			
			if (between(sstart, send, vstart) || between(sstart, send, vend) || between(vstart, vend, sstart)){
				return 0;
			} else {
				return sstart.compareTo(vstart);
			}
			
			
		} else {
			return chrcomp;
		}
	}
	
	public static int snvIntersect(MAFData s, MAFData v) {
		// the subject is an snv so we don't expand and require a base match
		int chrcomp = s.get("Chromosome").compareTo(v.get("Chromosome"));
		if (chrcomp == 0){
			Integer sstart = Integer.parseInt(s.get("Start_position")) - 5;
			Integer vstart = Integer.parseInt(v.get("Start_position"));
			int pcomp = sstart.compareTo(vstart);
			if (pcomp == 0){
				return s.get("Tumor_Sequence_Allele2").compareTo(v.get("Tumor_Sequence_Allele2"));
			} else {
				return pcomp;
			}
		} else {
			return chrcomp;
		}
	}

}
