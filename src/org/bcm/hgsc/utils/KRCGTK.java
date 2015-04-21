package org.bcm.hgsc.utils;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;

public class KRCGTK {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Options options = new Options();
		Parser parser = new BasicParser();
		options.addOption("a", true, "action to process");
	}

}
