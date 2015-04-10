package org.bcm.hgsc.cancer.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

public class VCFConverter {
	/**
     * @param args the command line arguments
     */
    private static void command() {
        System.out.println("Java -Xmx4g -jar vcfs.jar");
        System.out.println("Java -Xmx4g -jar vcfs.jar [-project <project-path> -output <output-path>");
        System.out.println("     project -- specify a project path");
        System.out.println("     output  -- specify a output path");
        System.out.println("example: Java -Xmx4g -jar vcfs.jar -proj /stornext/snfs6/cancer-analysis/thca/we -out vcf1");
    }

    public static void main(String[] args) throws ParseException {
    	
    	Options options = new Options();
		Parser parser = new BasicParser();
		
		
		HelpFormatter formatter = new HelpFormatter();
		CommandLine line = parser.parse(options, args);
		
		String[] reqargs = line.getArgs();
		
		if (reqargs.length < 2) {
			System.out.println("Required args not supplied.");
			formatter.printHelp("VCFConverter.jar projectDir mafs", options);
			System.exit(10);
		}
		String projdir = reqargs[0];
        String[] mafs = Arrays.copyOfRange(reqargs, 1, reqargs.length);
        convert(projdir, mafs);
        //convert(designPath, ".removed", outPath);
    }

    public static String getCanonicalPath(File file) {
        try {
            return file.getCanonicalPath();
        } catch (Exception e) {
        }
        return "";
    }

    public static void convert(String projdir, String[] mafs) {//as we or wgs' path

    	// first build a list of sites to remove from the vcf files
    	ArrayList<String> removedVariants = new ArrayList<String>();
    	String line;
    	String[] items;
    	ArrayList<String> chrs;
    	ArrayList<Integer> ids;
    	int i, j, k, chr, start, ref, var, tid;
    	for (i = 0; i < mafs.length; i++) {
    		try {
    			System.err.println("Search filtered file:" + mafs[i]);
    			BufferedReader reader = new BufferedReader(new FileReader(mafs[i]));
    			line = reader.readLine();
    			items = line.split("\t");
    			chr = start = tid = ref = var = -1;
    			for (j = 0; j < items.length; j++) {
    				line = items[j].toLowerCase();
    				if ("chromosome".startsWith(line)) {
    					chr = j;
    				} else if ("start_position".startsWith(line)) {
    					start = j;
    				} else if (line.indexOf("tumor") >= 0 && line.indexOf("barcode") >= 0) {
    					tid = j;
    				} else if ("reference_allele".startsWith(line)) {
    					ref = j;
    				} else if ("tumor_seq_allele2".startsWith(line)) {
    					var = j;
    				}
    			}
    			while ((line = reader.readLine()) != null) {
    				items = line.split("\t");
    				//line = items[tid].trim().toUpperCase() + "\t" + items[chr].trim().toUpperCase() + "\t" + items[start].trim() + "\t" + items[ref].trim().toUpperCase() + "\t" + items[var].trim().toUpperCase();
    				line = items[tid].trim().toUpperCase() + "\t" + items[chr].trim().toUpperCase() + "\t" + items[start].trim();
    				if (!removedVariants.contains(line)) {
    					removedVariants.add(line);
    				}
    			}
    			reader.close();
    		} catch (Exception e) {
    		}
    	}
    	convertToVCF(removedVariants, projdir + File.separator + "atlas-snp" + File.separator + "filtered");
    	convertToVCF(removedVariants, projdir + File.separator + "atlas-indel" + File.separator + "filtered");
    	if (removedVariants.size() > 0) {
    		System.err.println("************Variants not detected " + removedVariants.size() + "*************");
    		for (i = 0; i < removedVariants.size(); i++) {
    			System.err.println(removedVariants.get(i) + "\n");
    		}
    	}
    	//convertToVCF(sampleIds, removedVariants, designPath, outpath);qq

    }

    public static void convertToVCF(ArrayList<String> removedVariants, String filterpath) {
        File file = new File(filterpath);
        if (file.exists()) {
            boolean indel;
            int i, j, n1, n2, chr, pos, filter, ref, var, info;
            String line, fn, item, id;
            String[] items;
            //Object[] obs;
            //ArrayList<String> chrs;
            //ArrayList<Integer> ids;
            File[] files = file.listFiles();
            for (i = 0; i < files.length; i++) {
                if (!files[i].getName().startsWith(".") && !files[i].isHidden() && (line = getCanonicalPath(files[i])).endsWith(".vcf")) {
                    try {
                        indel = (line.indexOf(".indel.") > 0);
                        j = line.indexOf(File.separator + "atlas-");
                        fn = System.getProperty("user.dir");
                        if (j >= 0) {
                            fn += line.substring(j);
                        }
                        j = fn.lastIndexOf(File.separator + "filtered");
                        if (j > 0) {
                            fn = fn.substring(0, j);
                        }
                        File f = new File(fn);
                        if (!f.exists()) {
                            f.mkdirs();
                        }
                        j = line.lastIndexOf(File.separator);
                        fn += line.substring(j++);
                        System.err.println("Convert VCF to " + fn + " from " + line);
                        line = line.substring(j);
                        j = line.indexOf(".");
                        id = line.substring(0, j).toUpperCase();
                        //id = sampleIds.indexOf(line);
                        FileWriter writer = new FileWriter(fn);
                        BufferedReader reader = new BufferedReader(new FileReader(files[i]));
                        chr = pos = filter = ref = var = info = -1;
                        while ((line = reader.readLine()) != null && line.startsWith("##")) {
                            writer.write(line + "\n");
                            if (line.startsWith("##FILTER=")) {
                                writer.write("##FILTER=<ID=POSTCA,Description=\"Fail Post Carnac (Germline, Gene count, normal variant, allele fraction)\">\n");
                            }
                        }
                        if (line.startsWith("#")) {
                            writer.write(line + "\n");
                            items = line.substring(1).split("\t");
                            for (j = 0; j < items.length; j++) {
                                line = items[j].toUpperCase();
                                if ("CHROMOSOME".startsWith(line)) {
                                    chr = j;
                                } else if ("POSITION".startsWith(line)) {
                                    pos = j;
                                } else if ("FILTER".startsWith(line)) {
                                    filter = j;
                                } else if ("REF".startsWith(line)) {
                                    ref = j;
                                } else if ("ALT".startsWith(line)) {
                                    var = j;
                                } else if ("INFO".startsWith(line)) {
                                    info = j;
                                }
                            }
                        }
                        n1 = n2 = 0;
                        while ((line = reader.readLine()) != null) {
                            items = line.split("\t");
                            item = items[filter].toUpperCase();
                            if (item.equals("PASS")) {
                                n1++;
                                //item = id + "\t" + items[chr].trim().toUpperCase() + "\t" + (indel ? Integer.parseInt(items[pos].trim())+1 : items[pos].trim())+ "\t" + items[ref].trim().toUpperCase() + "\t" + items[var].trim().toUpperCase();
                                item = id + "\t" + items[chr].trim().toUpperCase() + "\t" + (indel && items[info].indexOf("VT=DEL") >= 0 ? Integer.parseInt(items[pos].trim()) + 1 : items[pos].trim());
                                if (removedVariants.contains(item)) {
                                    items[filter] = "POSTCA";
                                    removedVariants.remove(item);
                                    n2++;
                                    items[filter] = "POSTCA";
                                    for (j = 0; j < items.length - 1; j++) {
                                        writer.write(items[j] + "\t");
                                    }
                                    writer.write(items[j] + "\n");
                                } else {
                                    writer.write(line + "\n");
                                }                               
                            } else {
                                writer.write(line + "\n");
                            }
                        }
                        System.err.println("Converted PASS to POSTCA " + n2 + " out of " + n1);
                        reader.close();
                        writer.close();

                    } catch (Exception e) {
                    }
                }
            }
        }
    }
}
