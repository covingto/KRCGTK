package org.bcm.hgsc.utils;

public class RangUtils {
	
	public static boolean between(double a, double b, double q){
		if (q <= b & q >= a) { return true; }
		return false;
	}
	
	public static boolean intersects(int s1, int e1, int s2, int e2){
		if (between(s1, e1, s2) || between(s1, e1, e2) || between(s2, e2, s1)){ return true; }
		else { return false; }
	}
}
