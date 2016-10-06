package org.bcm.hgsc.utils;

import java.io.File;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.logging.Logger;
import java.util.Set;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SamReaderFactory.Option;

public class BAMInterface {
	private final File samfilereader;
	private final String sampleName;
	private final String sampleType;
	private static final Logger log = Logger.getLogger(BAMInterface.class.getName());
	
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
	
	public SamReader getSamfilereader(final Set<Option> enabledOptions, final Set<Option> disabledOptions){
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
		for (Option o : enabledOptions){
			readerFactory.enable(o);
		}
		for (Option o : disabledOptions){
			readerFactory.disable(o);
		}
		for (Entry<String, String> entry : Settings.env.entrySet()){
			final String key = entry.getKey();
			final String value = entry.getValue();
			switch(key){
			case "SAM_VALIDATION_STRINGENCY":
				switch(value){
				case "LENIENT":
					readerFactory.validationStringency(ValidationStringency.LENIENT);
					break;
				case "SILENT":
					readerFactory.validationStringency(ValidationStringency.SILENT);
					break;
				case "STRICT":
					readerFactory.validationStringency(ValidationStringency.STRICT);
					break;
				case "NULL":
					break;
				default:
					entry.setValue("NULL");
					log.warning("SAM_VALIDATION_STRINGENCY is not an appropriate value");
					break;
				}
				break;
			default:
				break;
			}
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
