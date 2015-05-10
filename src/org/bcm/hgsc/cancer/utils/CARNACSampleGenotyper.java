package org.bcm.hgsc.cancer.utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.lang3.StringUtils;
import org.bcm.hgsc.utils.AlleleResolver.AlleleSet;
import org.bcm.hgsc.utils.BAMUtils.ConformedRead;
import org.bcm.hgsc.utils.ByteContainer;
import org.bcm.hgsc.utils.Utils;


/**
 * Primary class for genotyping (that is returning allele counts and info at genomic locations) for a given sample
 * @author covingto
 *
 */
public class CARNACSampleGenotyper extends SampleGenotyper{
	//private final String sampleName;
	//private final File bamFilePath;
	//private final SAMFileReader samReader;
	private static Logger log = Logger.getLogger(AlleleSet.class.getName());
	public static int ENDBUFFER = 15;
	public static int MAXMAPQ = 60;
	static {
		headerlines.add(new VCFHeaderLine("Genotyper", "org.bcm.hgsc.cancer.CARNACSampleGenotyper"));
		headerlines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		headerlines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.String, "Read Depth"));
		headerlines.add(new VCFInfoHeaderLine("OC", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Alleles called by each caller source listed as the source file of the initial call.  These may be modified during allele resolution so the origingal position and allele for each call is listed."));
		headerlines.add(new VCFFormatHeaderLine("AC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total coverage for each allele"));
		headerlines.add(new VCFFormatHeaderLine("AQC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total coverage at Q20 for each allele"));
		headerlines.add(new VCFFormatHeaderLine("ASQ", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Sum of the qualities for each allele"));
		headerlines.add(new VCFFormatHeaderLine("FC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total forward for each allele"));
		headerlines.add(new VCFFormatHeaderLine("RC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total reverse for each allele"));
		headerlines.add(new VCFFormatHeaderLine("FQC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total forward coverage at Q20 for each allele"));
		headerlines.add(new VCFFormatHeaderLine("RQC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total reverse coverage at Q20 for each allele"));
		headerlines.add(new VCFFormatHeaderLine("MR", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total mid-read coverage for each allele"));
		headerlines.add(new VCFFormatHeaderLine("MMQ", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Maximum mapping quality of reads supporting each allele"));
		headerlines.add(new VCFFormatHeaderLine("MMQC", 1, VCFHeaderLineType.Integer, "Number of reads at this position that achieved maximum mapping quality"));
		headerlines.add(new VCFFormatHeaderLine("MRQ", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total mid-read Q20 allele count"));
		headerlines.add(new VCFFormatHeaderLine("R1C", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total read1 for each allele"));
		headerlines.add(new VCFFormatHeaderLine("R2C", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total read2 for each allele"));
		headerlines.add(new VCFFormatHeaderLine("R1QC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total read1 coverage at Q20 for each allele"));
		headerlines.add(new VCFFormatHeaderLine("R2QC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Total read2 coverage at Q20 for each allele"));
		
	}

	/**
	 * Genotyper for a sample at a position.
	 * For the set of alleles input each conformed read is tested against the allele to see if there is a match.
	 * Conformed reads are reduced to remove insertion placeholders prior to matchup.  Deletion placeholders are retained.
	 * Additionally local "realignment" occurs so that the alleles tested need not be in phase with one another to trigger a match.
	 * Graphically;
	 * Ex 1;
	 * Reference:		[	A	C	-	-	]
	 * Allele : 		[	A	T	C	G 	]
	 * Conformed Read:	[	A	T	C	G	]
	 * Conformed Reduc:	[	A	T	C	G	]
	 * Result: match
	 * Ex 2;
	 * Reference:		[	A	T	C	G	]
	 * Allele:			[	A	A	]
	 * Conformed Read:	[	A	-	-	A	]
	 * Conformed Reduc:	[	A	A	]
	 * Result: match
	 * Ex 3;
	 * Reference:		[	A	T	C	G	]
	 * Allele:			[	A	G	G	]
	 * Conformed Read1:	[	A	G	-	G	]
	 * Conformed Redu1:	[	A	G	G	]
	 * Conformed Read2:	[	A	-	G	G	]
	 * Conformed Redu2:	[	A	G	G	]
	 * Result 1: match
	 * Result 2: match
	 * Example 4;
	 * Reference:		[	C	G	G	G	G	T	]
	 * Allele:			[	C	G	G	G	A	]
	 * Conformed Read1:	[	C	-	G	G	G	A	]
	 * Conformed Redu1:	[	C	G	G	G	A	]
	 * Conformed Read2:	[	C	G	G	G	A	-	]
	 * Conformed Redu2:	[	C	G	G	G	A	]
	 * Result 1: match
	 * Result 2: match
	 * Note: Allele simplification will later reduce the allele to GGT>GA for both cases with an increment in position from i to i+3
	 * 
	 * 
	 * @param alleles
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	@Override
	public Genotype genotype(String sampleName, AlleleSet alleleset, List<ConformedRead> reads){
		// get some constants for this function
		final int sliceStart = alleleset.getSliceStart();
		final int sliceEnd = alleleset.getSliceEnd();
		final int start = alleleset.getStart();
		final int end = alleleset.getEnd();
		final int leftSliceOffset = start - sliceStart;
		final int rightSliceOffset = sliceEnd - end;
		final List<Allele> alleles = new ArrayList<Allele>(alleleset.getAlleles());

		// we will be filling in data for each allele for the following values
		final int arraySize = alleles.size();
		int[] forwardCoverage = new int[arraySize];
		int[] reverseCoverage = new int[arraySize];
		int[] forwardQ20Coverage = new int[arraySize];
		int[] reverseQ20Coverage = new int[arraySize];
		int[] midReadCoverage = new int[arraySize];
		int[] alleleCoverage = new int[arraySize];
		int[] alleleQ20Coverage = new int[arraySize];
		int[] alleleSumQual = new int[arraySize];
		int[] alleleMaxMapQual = new int[arraySize];
		int[] midReadQ20Coverage = new int[arraySize];
		int[] r1Coverage		=	new int[arraySize];
		int[] r2Coverage		=	new int[arraySize];
		int[] r1Q20Coverage		=	new int[arraySize];
		int[] r2Q20Coverage		=	new int[arraySize];
		int totalCoverage = 0;
		int maxMapQualCount = 0;

		// genotypeDP and genotypeMaxMapQualCount is determined from each conformed read
		for (ConformedRead cr : reads){
			// if the conformed read does not actually contain the allele, then we continue, this will reduce the actual number of alleles that we can recover
			if (cr.readStart() > sliceStart || cr.readEnd() < sliceEnd){ 
				continue; 
			}
			totalCoverage++;
			if (cr.getMapQuality() >= MAXMAPQ){
				maxMapQualCount++;
			}
			ByteContainer seqAllele = cr.getSeqAllele(sliceStart, sliceEnd);
			ByteContainer seqQuality = cr.getSeqQuality(sliceStart, sliceEnd);
			if (seqAllele.bytes.length - rightSliceOffset < leftSliceOffset){
				// this can apparently rarely happen
				log.log(Level.FINE, "Error processing allele for a read: Sequence: " + new String(seqAllele.bytes) + " leftSliceOffset: " + leftSliceOffset + " rightSliceOffset: " + rightSliceOffset);
				continue;
			}
			for (int i = 0; i < arraySize; i++){
				final Allele thisAllele = alleles.get(i);
				if (thisAllele.basesMatch(Arrays.copyOfRange(seqAllele.bytes, leftSliceOffset, seqAllele.bytes.length - rightSliceOffset))){
					alleleCoverage[i]++;
					// TODO: for quality we would like to indicate the quality of the first non-reference base
					// for now report the minimum quality of the allele
					int qual = seqQuality.bytes[0];
					for (int j = 0; j < seqQuality.bytes.length; j++){
						qual = Math.min(qual, seqQuality.bytes[j]);
					}
					final int deltaStart = start - cr.readStart();
					final int deltaEnd = cr.readEnd() - end;
					final int minEndDist = Math.min(deltaStart, deltaEnd);
					alleleSumQual[i] += qual;
					if (qual > 20){
						alleleQ20Coverage[i]++;
					}
					if (minEndDist > ENDBUFFER ){
						midReadCoverage[i]++;
						if (qual > 20){
							midReadQ20Coverage[i]++;
						}
					}
					if (cr.isForward()){
						forwardCoverage[i]++;
						if (qual > 20){
							forwardQ20Coverage[i]++;
						}
					} else {
						reverseCoverage[i]++;
						if (qual > 20){
							reverseQ20Coverage[i]++;
						}
					}
					if (cr.getMapQuality() > alleleMaxMapQual[i]){
						alleleMaxMapQual[i] = cr.getMapQuality();
					}
					if (cr.isRead1()){
						r1Coverage[i]++;
						if (qual > 20){
							r1Q20Coverage[i]++;
						}
					} else {
						r2Coverage[i]++;
						if (qual > 20){
							r2Q20Coverage[i]++;
						}
					}
					break;
				}
			}
		}

		GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName, alleles);
		genotypeBuilder.DP(totalCoverage);
		genotypeBuilder.attribute("AC",  StringUtils.join(Utils.intArrayToIntegerList(alleleCoverage), ","));
		genotypeBuilder.attribute("AQC", StringUtils.join(Utils.intArrayToIntegerList(alleleQ20Coverage), ","));
		genotypeBuilder.attribute("ASQ", StringUtils.join(Utils.intArrayToIntegerList(alleleSumQual), ","));
		genotypeBuilder.attribute("FC",  StringUtils.join(Utils.intArrayToIntegerList(forwardCoverage), ","));
		genotypeBuilder.attribute("RC",  StringUtils.join(Utils.intArrayToIntegerList(reverseCoverage), ","));
		genotypeBuilder.attribute("FQC", StringUtils.join(Utils.intArrayToIntegerList(forwardQ20Coverage), ","));
		genotypeBuilder.attribute("RQC", StringUtils.join(Utils.intArrayToIntegerList(reverseQ20Coverage), ","));
		genotypeBuilder.attribute("MR",  StringUtils.join(Utils.intArrayToIntegerList(midReadCoverage), ","));
		genotypeBuilder.attribute("MMQ", StringUtils.join(Utils.intArrayToIntegerList(alleleMaxMapQual), ","));
		genotypeBuilder.attribute("MMQC", maxMapQualCount);
		genotypeBuilder.attribute("MRQ", StringUtils.join(Utils.intArrayToIntegerList(midReadQ20Coverage), ","));
		genotypeBuilder.attribute("R1C", StringUtils.join(Utils.intArrayToIntegerList(r1Coverage), ","));
		genotypeBuilder.attribute("R2C", StringUtils.join(Utils.intArrayToIntegerList(r2Coverage), ","));
		genotypeBuilder.attribute("R1QC", StringUtils.join(Utils.intArrayToIntegerList(r1Q20Coverage), ","));
		genotypeBuilder.attribute("R2QC", StringUtils.join(Utils.intArrayToIntegerList(r2Q20Coverage), ","));
		
		return genotypeBuilder.make();
	}
}
