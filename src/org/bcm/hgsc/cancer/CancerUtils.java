package org.bcm.hgsc.cancer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.LinkedList;


public class CancerUtils {
	static String[][] readPairsFile(File pairs) throws IOException{
		BufferedReader reader = new BufferedReader(new FileReader(pairs));
		LinkedList<String[]> rows = new LinkedList<String[]>();
		String line;
		while ((line = reader.readLine()) != null){
			rows.add(line.split("\t"));
		}
		reader.close();
		return (String[][]) rows.toArray();
	}
	
	public class DirectoryFilter implements FilenameFilter{
		private String dirname;
		
		public DirectoryFilter(String dirname){
			super();
			this.dirname = dirname;
		}

		@Override
		public boolean accept(File f, String s) {
			if (s == this.dirname){
				return true;
			}
			return false;
		}
	}
	
//	public File[] findFile(File directory, String pattern){
//		public static class Finder
//			extends SimpleFileVisitor<Path> {
//
//			private final PathMatcher matcher;
//			private int numMatches = 0;
//
//			Finder(String pattern) {
//				matcher = FileSystems.getDefault()
//						.getPathMatcher("glob:" + pattern);
//			}
//
//			// Compares the glob pattern against
//			// the file or directory name.
//			void find(Path file) {
//				Path name = file.getFileName();
//				if (name != null && matcher.matches(name)) {
//					numMatches++;
//					System.out.println(file);
//				}
//			}
//
//			// Prints the total number of
//			// matches to standard out.
//			void done() {
//				System.out.println("Matched: "
//						+ numMatches);
//			}
//
//			// Invoke the pattern matching
//			// method on each file.
//			@Override
//			public FileVisitResult visitFile(Path file,
//					BasicFileAttributes attrs) {
//				find(file);
//				return CONTINUE;
//			}
//
//			// Invoke the pattern matching
//			// method on each directory.
//			@Override
//			public FileVisitResult preVisitDirectory(Path dir,
//					BasicFileAttributes attrs) {
//				find(dir);
//				return CONTINUE;
//			}
//
//			@Override
//			public FileVisitResult visitFileFailed(Path file,
//					IOException exc) {
//				System.err.println(exc);
//				return CONTINUE;
//			}
//		}
//	}
}
