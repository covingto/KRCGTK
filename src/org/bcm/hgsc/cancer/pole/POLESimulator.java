package org.bcm.hgsc.cancer.pole;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.bcm.hgsc.cancer.utils.AlternativeHypothesis;
import org.bcm.hgsc.cancer.utils.BinomialTest;

public class POLESimulator {
	// public BinomialTest bt = new BinomialTest();
	public enum Strand{
		POS, NEG;
	}
	public enum Type{
		TCT, AGA, WT;
	}
	public enum Rep{
		S, US;
	}

	// =============================
	public class POLE{
		private int mPos;
		// will run simulation on 100 cell divisions so even through the rate is 1600 / 1000000 we will do this 100 times
		private float mMutrate = (float) 10 / (float) 1000000;
		private int mStartTime;
		private Strand mStrand;

		// =============================
		public POLE(int t, int p, Strand s){
			this.mPos = p;
			this.mStartTime = t;
			this.mStrand = s;
		}

		// =============================
		public boolean advance(int t){
			// check to see if the POLE has started yet
			if (t < mStartTime){
				return true;
			}
			if (mStartTime >= POLESimulator.mGenomeSize){
				return false;
			}
			// advance the position
			if (mStrand.equals(Strand.POS)){
				mPos++;
			} else {
				mPos--;
			}
			// check to see if we are still in the region
			if (mPos >= 0 && mPos < rep.length){
				// check to see if the site has been synthesized yet
				if (POLESimulator.rep[mPos].equals(Rep.S)){
					//System.out.println("Base replicated stopping.");
					return false;
				} else if (!POLESimulator.type[mPos].equals(Type.WT)){
					return true; // already mutated in an earlier iteration.
				} else {
					POLESimulator.rep[mPos] = Rep.S;
					// fire and update the mutation
					if (POLESimulator.mRand.nextFloat() < mMutrate){
						if (mStrand.equals(Strand.POS)){
							if (mPos % 2 == 0){ 
								//System.out.println("Mutating " + mPos + " to tct");
								POLESimulator.type[mPos] = Type.TCT; 
							}
						} else {
							if (mPos % 2 == 1){ 
								//System.out.println("Mutating " + mPos + " to aga");
								POLESimulator.type[mPos] = Type.AGA; 
							}
						}
					}
				}
			} else {
				return false;
			}
			return true;

		}
	}

	// =============================
	public static int mSimulationNumber = 0;
	public static int mWindowSize = 20000;
	public static int mJumpSize = 1000;
	public static int	mGenomeSize = 25000000;
	public static Rep[] rep = 	new Rep[mGenomeSize];
	public static Type[] type =	new Type[mGenomeSize];
	public static Random mRand = new Random();


	// =============================
	public static void resetGenome(){
		System.out.println("Returning to wildtype genome");
		Arrays.fill(rep, Rep.US);
		Arrays.fill(type, Type.WT);
	}

	// =============================
	public void runRandomORI(int numORI, int oriStartRandMax){	
		Arrays.fill(rep, Rep.US);
		POLE[] lPOLEArray = new POLE[numORI * 2];
		// set up the POLE's
		for (int i = 0; i < lPOLEArray.length; i++){
			final int l1P = mRand.nextInt(mGenomeSize);
			final int l1StartT = mRand.nextInt(oriStartRandMax);
			lPOLEArray[i] = new POLE(l1StartT, l1P, Strand.POS);
			i++;
			lPOLEArray[i] = new POLE(l1StartT, l1P, Strand.NEG);
		}
		runSimulation(lPOLEArray);
	}

	// ============================
	public void runSpacedTimedORI(int[] oriPos, int[] oriTime) throws Exception{
		Arrays.fill(rep, Rep.US);
		POLE[] lPOLEArray = new POLE[oriPos.length * 2];
		if (oriPos.length != oriTime.length){
			throw new Exception("oriPos length != oriTime length");
		}
		for (int i = 0; i < oriPos.length; i++){
			final int l1P = oriPos[i];
			final int l1StartT = oriTime[i];
			lPOLEArray[i*2] = new POLE(l1StartT, l1P, Strand.POS);
			lPOLEArray[i*2 + 1] = new POLE(l1StartT, l1P, Strand.NEG);
		}
		runSimulation(lPOLEArray);
	}


	// ============================
	public static void runSimulation(POLE[] lPOLEArray){
		System.out.println("Running replication sumulation " + mSimulationNumber);
		mSimulationNumber++;
		int t = 0; // t represents the time point in the experiment

		while ( t < mGenomeSize ){
			// boolean[] toRemove = new boolean[lPOLEArray.length];
			// if (t % 1000000 == 0){ System.out.println("  Position " + t); }
			int doneCount = 0;
			for (int j = 0; j < lPOLEArray.length; j++){
				final POLE l2POLE = lPOLEArray[j];
				if (!l2POLE.advance(t)){ doneCount++;}
			}
			if (doneCount == lPOLEArray.length){
				// they are all done so we can end class early
				System.out.println("Exiting replication after " + t + " iterations.");
				break;
			}
			t++;
		}
		System.out.println("Replication done");
	}

	// ============================
	public static List<Double> binomWindow(){
		List<Double> pvals = new ArrayList<Double>(mGenomeSize / mJumpSize);
		for (int i = 0; i < (mGenomeSize - mWindowSize); i = i + mJumpSize){
			int tctCount = 0;
			int agaCount = 0;
			//System.out.println("Types");
			for (int j = i; j < i + mWindowSize; j++){
				//System.out.print(type[j] + ", ");
				if (type[j].equals(Type.TCT)){ 
					tctCount++;
				} else if (type[j].equals(Type.AGA)){
					agaCount++;
				}
			}
			// System.out.println();
			int tctDom = tctCount > agaCount ? 1 : -1;
			try{
				double btr = Math.min(1, BinomialTest.binomialTest(tctCount + agaCount, tctCount, 0.5, AlternativeHypothesis.TWO_SIDED));
				//System.out.println(btr + "\t" + tctDom + "\t" + agaCount + "\t" + tctCount);
				pvals.add(tctDom * btr);
			} catch (Exception e){
				pvals.add(new Double(tctDom));
			}
			//System.out.println("Found high tct");
		}
		return pvals;

	}

	// ============================
	public static Double binomP(boolean[] tctvaga){
		int tctCount = 0;
		for (int i = 0; i < tctvaga.length; i++){
			if (tctvaga[i]){ tctCount++; }
		}
		return BinomialTest.binomialTest(tctvaga.length, tctCount, 0.5, AlternativeHypothesis.TWO_SIDED);
	}

	// ============================
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		// simulate random data 100X
		POLESimulator sim = new POLESimulator();
		int[] startTime = new int[100];
		int[] startPos = new int[100];
		Arrays.fill(startTime, 0);
		int[] startTimeA = Arrays.copyOf(startTime, startTime.length);
		int[] startTimeB = Arrays.copyOf(startTime, startTime.length);
		// simulate equal spacing all same timing
		// placement 1000 equally spaced ori's across the genome
		for (int i = 0; i < 100; i++){
			int genomeMultiplier = mGenomeSize / 100;
			startPos[i] = i * genomeMultiplier;
		}


		for (int i = 0; i < startTimeB.length; i++){
			if (i % 4 == 0){
				startTimeB[i] = mGenomeSize; // never starts
			}
		}
		if (args[0].equals("rand")){
			File randomOutputFile = new File("randomPlaceTime.txt");
			BufferedWriter randomWriter = new BufferedWriter(new FileWriter(randomOutputFile));

			for (int s1 = 0; s1 < 10; s1++){
				POLESimulator.resetGenome();
				for (int s1r = 0; s1r < 100; s1r++){
					sim.runRandomORI(100, 10000);
				}
				List<Double> bwr = POLESimulator.binomWindow();
				int ln = 1;
				System.out.println("Writing trial results");
				for (Double bt : bwr){
					randomWriter.write(s1 + "\t" + ln + "\t" + bt.toString());
					randomWriter.newLine();
					ln++;
				}
			}
			randomWriter.flush();
			randomWriter.close();
		} else if (args[0].equals("clampPT")){

			File equalOutputFile = new File("equalPlaceTime.txt");
			BufferedWriter equalWriter = new BufferedWriter(new FileWriter(equalOutputFile));
			for (int s1 = 0; s1 < 10; s1++){
				POLESimulator.resetGenome();
				for (int s1r = 0; s1r < 100; s1r++){
					sim.runSpacedTimedORI(startPos, startTime);
				}
				List<Double> bwr = POLESimulator.binomWindow();
				int ln = 1;
				System.out.println("Writing trial results");
				for (Double bt : bwr){
					equalWriter.write(s1 + "\t" + ln + "\t" + bt);
					equalWriter.newLine();
					ln++;
				}
			}
			equalWriter.flush();
			equalWriter.close();
		} else if (args[0].equals("clampPsrandT")){
			// simulate equal spacing all same timing
			// placement 1000 equally spaced ori's across the genome

			File equalPSwitchTOutputFile = new File("equalPlaceSwitchTime.txt");
			BufferedWriter equalPSwitchTWriter = new BufferedWriter(new FileWriter(equalPSwitchTOutputFile));
			for (int s1 = 0; s1 < 10; s1++){
				POLESimulator.resetGenome();
				for (int s1r = 0; s1r < 100; s1r++){
					sim.runSpacedTimedORI(startPos, s1 % 2 == 0 ? startTimeA : startTimeB);
				}
				List<Double> bwr = POLESimulator.binomWindow();
				int ln = 1;
				System.out.println("Writing trial results");
				for (Double bt : bwr){
					equalPSwitchTWriter.write(s1 + "\t" + ln + "\t" + bt);
					equalPSwitchTWriter.newLine();
					ln++;
				}
			}
			equalPSwitchTWriter.flush();
			equalPSwitchTWriter.close();
		}
	}

}
