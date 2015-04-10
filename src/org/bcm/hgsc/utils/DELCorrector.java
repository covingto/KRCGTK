package org.bcm.hgsc.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DELCorrector {
	final static byte n = 'N';
	private static byte lowQual = 33;
	public static void main(String[] args) throws CloneNotSupportedException, IOException {
		if (args.length < 4){
			System.out.println("Useage: DELCorector.jar : inputBAM.bam outputBAM.bam correctorStart correctorEnd");
			System.out.println("Used to repair deletion error caused by machine errors in data aquisition.  correctorStart and correctorEnd reflect the boundries in reads within which D calls are replaced with N, quality set to int 33 (phred 0)");
			System.exit(10);
		}
		String infile = args[0];
		String outfile = args[1];
		Integer startRange = Integer.parseInt(args[2]);
		Integer endRange = Integer.parseInt(args[3]);
		
		final SamReader infileReader = SamReaderFactory.makeDefault().open(new File(infile));
		final SAMFileWriter outfileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(infileReader.getFileHeader(),
				true, new File(outfile));
		
		for (final SAMRecord rec : infileReader){
			SAMRecord fixedRec = repair(rec, startRange, endRange);
			outfileWriter.addAlignment(fixedRec);
		}
		
		infileReader.close();
		outfileWriter.close();

	}
	
	public static SAMRecord repair(SAMRecord rec, int startRange, int endRange) throws CloneNotSupportedException{
		final byte[] orgSeq = rec.getReadBases();
		final byte[] orgQual = rec.getBaseQualities();
		final byte[] orgOrgQual = rec.getOriginalBaseQualities();
		final Cigar orgCigar = rec.getCigar();
		
		int readDist = 0;
		int offendingCigElement = 0; // place holder to find offending cigar elements
		int numToFill = 0;
		// isolate the position of the del
		for (final CigarElement cigEl : orgCigar.getCigarElements()){
			final int cigElLen = cigEl.getLength();
			final CigarOperator cigElOp = cigEl.getOperator();
			if (cigElOp == CigarOperator.MATCH_OR_MISMATCH || cigElOp == CigarOperator.SOFT_CLIP || cigElOp == CigarOperator.INSERTION){
				readDist += cigElLen;
				offendingCigElement++;
			} else if (readDist > startRange && readDist < endRange && cigElOp == CigarOperator.DELETION){ // we got one!!!
				
				break;
			}
		}
		if (readDist >= endRange && offendingCigElement > 0 && offendingCigElement < (orgCigar.numCigarElements() - 1)){ return rec; } // we don't need to fix anything
		SAMRecord fixedRec = (SAMRecord) rec.clone();
		Cigar fixedCigar = new Cigar(resolveCigElements(orgCigar, offendingCigElement));
		byte[] fixedSeq = resolveSeq(readDist, numToFill, orgSeq);
		byte[] fixedQual = resolveQual(readDist, numToFill, orgQual);
		byte[] fixedOrgQual = resolveQual(readDist, numToFill, orgOrgQual);
		fixedRec.setAttribute("MD", null);
		fixedRec.setBaseQualities(fixedQual);
		fixedRec.setReadBases(fixedSeq);
		fixedRec.setOriginalBaseQualities(fixedOrgQual);
		fixedRec.setCigar(fixedCigar);
		return fixedRec;
	}
	
	public static byte[] resolveSeq(int readDist, int numToFill, byte[] orgSeq){
		byte[] newSeq = new byte[orgSeq.length + numToFill];
		int i = 0;
		for (i = 0; i < readDist; i++){
			newSeq[i] = orgSeq[i];
		}
		for (i = i; i < readDist + numToFill; i++){
			newSeq[i] = n;
		}
		for (i = i; i < orgSeq.length + numToFill; i++){
			newSeq[i] = orgSeq[i-numToFill];
		}
		return newSeq;
	}
	
	public static byte[] resolveQual(int readDist, int numToFill, byte[] orgQual){
		byte[] newQual = new byte[orgQual.length + numToFill];
		int i = 0;
		for (i = 0; i < readDist; i++){
			newQual[i] = orgQual[i];
		}
		for (i = i; i < readDist + numToFill; i++){
			newQual[i] = lowQual ;
		}
		for (i = i; i < orgQual.length + numToFill; i++){
			newQual[i] = orgQual[i-numToFill];
		}
		return newQual;
	}
	
	public static List<CigarElement> resolveCigElements(Cigar orgCigar, int offending){
		List<CigarElement> newElements = new ArrayList<CigarElement>();
		for (int i = 0; i < offending - 1; i++){
			final CigarElement cigEl = orgCigar.getCigarElement(i);
			newElements.add(cigEl);
			
		}
		// make the fix
		int ciglen = 
				orgCigar.getCigarElement(offending - 1).getLength() + 
				orgCigar.getCigarElement(offending - 1).getLength() + 
				orgCigar.getCigarElement(offending - 1).getLength();
		CigarElement fixedEl = new CigarElement(ciglen, CigarOperator.MATCH_OR_MISMATCH);
		newElements.add(fixedEl);
		// fill in the end
		for (int i = offending + 2; i < orgCigar.numCigarElements(); i++){
			final CigarElement cigEl = orgCigar.getCigarElement(i);
			newElements.add(cigEl);
		}
		return newElements;
	}
}
