package org.bcm.hgsc.cancer.sv;

public interface JunctionSourceReader {

	public boolean hasNext();

	public JRecord next();
	
	public void close();

}
