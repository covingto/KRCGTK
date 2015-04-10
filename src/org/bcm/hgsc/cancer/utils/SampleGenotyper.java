package org.bcm.hgsc.cancer.utils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.bcm.hgsc.utils.AlleleResolver.AlleleSet;
import org.bcm.hgsc.utils.BAMUtils.ConformedRead;

public abstract class SampleGenotyper {
	protected static Set<VCFHeaderLine> headerlines = new HashSet<VCFHeaderLine>();
	public static Set<VCFHeaderLine> getHeaderLines(List<String> samples, File sampleInfo){
		Set<VCFHeaderLine> newlines = new HashSet<VCFHeaderLine>();
		newlines.addAll(headerlines);
		if (sampleInfo != null && sampleInfo.canRead()){
			// TODO: read the sample info file and fill in the meta tags
		}
		return headerlines;
	}
	public static Set<VCFHeaderLine> getHeaderLines(){
		return headerlines;
	}
	public abstract Genotype genotype(String sampleName, AlleleSet alleleset, List<ConformedRead> reads);
}
