package org.bcm.hgsc.cancer.utils;


public abstract class WalkAcceptor {
	public abstract void add(String chr, long posa, long posb, double value) throws Exception;
}
