package org.bcm.hgsc.cancer.sv;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class NonRefSignal {
	private class NonRefSignalContainer {
		private int obs = 0;
		private List<Long> signalMM = new ArrayList<Long>();
		private List<Long> signalIM = new ArrayList<Long>();
		private List<Long> positions = new ArrayList<Long>();
		
		public void add(long pos, long mm, long im){
			// add to the signal and positions lists
			signalMM.add(mm);
			signalIM.add(im);
			positions.add(pos);
			
			// calculate the new average
			obs += 1;
		}
		
		private double average(List<Long> signal){
			long a = 0;
			for (Long v : signal){
				a += v;
			}
			return a / (double) obs;
		}
		
		public double averageMM(){
			return average(this.signalMM);
		}
		
		public double averageIM(){
			return average(this.signalIM);
		}
	}
	
	private Map<String, NonRefSignalContainer> containers = new HashMap<String, NonRefSignalContainer>();
	private final String[] chrs;
	
	public NonRefSignal(List<String> chromosomes){
		chrs = new String[chromosomes.size()];
		for (int i = 0; i < chromosomes.size(); i++){
			final String s = chromosomes.get(i);
			chrs[i] = s;
			containers.put(s, new NonRefSignalContainer());
		}
	}
	
	public void add(String chr, long posa, long posb, double signalMM, double signalIM){
		
	}
}
