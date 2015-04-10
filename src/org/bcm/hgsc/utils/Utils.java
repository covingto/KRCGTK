package org.bcm.hgsc.utils;

import java.io.File;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Utils {
	private static Logger log = Logger.getLogger(Utils.class.getName());
	public static File[] filesFromStrings(String[] filepaths){
		File[] newfiles = new File[filepaths.length];
		for (int i = 0; i < filepaths.length; i++){
			newfiles[i] = new File(filepaths[i]);
		}
		return newfiles;
	}
	public static String reverseComplement(String s){
		String rev = reverse(s);
		return complement(rev);
	}

	public static String reverse(String s){
		String r = "";
		for (int k = 0; k < s.length(); k++){
			r = s.charAt(k) + r;
		}
		return r;
	}

	public static String complement(String s){
		String r = "";
		for (int k = 0; k < s.length(); k++){
			if (s.charAt(k) == 'A'){
				r = r + 'T';
			} else if (s.charAt(k) == 'T'){
				r = r + 'A';
			} else if (s.charAt(k) == 'C'){
				r = r + 'G';
			} else if (s.charAt(k) == 'G'){
				r = r + 'C';
			} else if (s.charAt(k) == 'a'){
				r = r + 't';
			} else if (s.charAt(k) == 't'){
				r = r + 'a';
			} else if (s.charAt(k) == 'c'){
				r = r + 'g';
			} else if (s.charAt(k) == 'g'){
				r = r + 'c';
			} else {
				r = r + 'N';
			}
		}
		return r;
	}

	public static byte[] bytesToBytes(Byte[] bytes){
		byte[] result = new byte[bytes.length];
		for (int i = 0; i < bytes.length; i++){
			result[i] = bytes[i];
		}
		return result;
	}

	public static List<Integer> intArrayToIntegerList(int[] ints){
		List<Integer> newints = new ArrayList<Integer>(ints.length);
		for (int i = 0; i < ints.length; i++){
			newints.add(ints[i]);
		}
		return newints;
	}

	public static List<Byte> byteArrayToByteList(byte[] bytes){
		List<Byte> newints = new ArrayList<Byte>(bytes.length);
		for (int i = 0; i < bytes.length; i++){
			newints.add(bytes[i]);
		}
		return newints;
	}

	public static <T> Map<T, Integer> countOccurrences(Collection<T> list){
		Map<T, Integer> occurrenceMap = new HashMap<T, Integer>();

		for (T obj: list){
			Integer numOccurrence = occurrenceMap.get(obj);
			if(numOccurrence == null){
				//first count
				occurrenceMap.put(obj, 1);
			} else{
				occurrenceMap.put(obj, numOccurrence+ 1);
			}
		}

		return occurrenceMap;
	}
	public static void jmap(File file) {
		String name = ManagementFactory.getRuntimeMXBean().getName();
		String[] str = name.split("@");
		try {
			Runtime.getRuntime().exec("jmap -dump:file=" + file + " " + str[0]);
		} catch (IOException e) {
			log.log(Level.WARNING, "Exception running jmap", e);
		}
	}
}
