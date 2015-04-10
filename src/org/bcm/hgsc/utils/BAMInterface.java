package org.bcm.hgsc.utils;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReaderFactory.Option;

import java.io.File;
import java.util.EnumSet;

public class BAMInterface {
	private final File samfilereader;
	private final String sampleName;
	private final String sampleType;
	
	public BAMInterface(File samfilereader, String sampleName, String sampleType) throws Exception{
		this.samfilereader = samfilereader;
		this.sampleName = sampleName;
		this.sampleType = sampleType;
		if (this.samfilereader == null){
			throw new Exception("SAM file reader is null");
		}
	}
	
	public BAMInterface(String samfile, String sampleName, String sampleType) throws Exception{
		this(new File(samfile), sampleName, sampleType);
	}

	/**
	 * Generates a {@link SamReader} for the file backed by this interface.  The {@link SamReader} is opened but not closed by the interface.
	 * Note that by default this implementation ENABLES {@link SamReaderFactory.Option#DONT_MEMORY_MAP_INDEX} and DISABLES {@link SamReaderFactory.Option#CACHE_FILE_BASED_INDEXES}.
	 * Use an alternate constructor to specify the enabled and disabled options for this reader.
	 * @return
	 */
	public SamReader getSamfilereader() {
		return getSamfilereader(EnumSet.of(Option.DONT_MEMORY_MAP_INDEX), EnumSet.of(Option.CACHE_FILE_BASED_INDEXES));
	}
	
	public SamReader getSamfilereader(final EnumSet<Option> enabledOptions, final EnumSet<Option> disabledOptions){
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
		for (Option o : enabledOptions){
			readerFactory.enable(o);
		}
		for (Option o : disabledOptions){
			readerFactory.disable(o);
		}
		return readerFactory.open(samfilereader);
	}

	public String getSampleName() {
		return sampleName;
	}

	public String getSampleType() {
		return sampleType;
	}
	
}
