package org.bcm.hgsc.utils;

public class NonRefBaseCounter {
	/**
	 * Idea is to take a bam file and count the number of non-reference bases.  This counting is done in a "mapreduce" manner whereby
	 * the program runs in multiple threads, each processing reads from a block of reads and counting the number of of non-reference
	 * positions.  These are output in the form A>C, A>G, ... a>c, a>g, ..., ins, del
	 * The ultimate output of this program is a file with two lines of text, a header line and a data line.
	 * @param args
	 */
	public static void main(String[] args) {
		
	}

}
