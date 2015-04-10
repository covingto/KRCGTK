package org.bcm.hgsc.cancer.pacbio;

public class BlankSyntheticReference extends SyntheticReference {
	public BlankSyntheticReference(){
		super();
		this.setCantClip();
	}
	
	public SyntheticReference clip(Integer leftclip) throws ExcessiveClipException, DoubleClipException {
		return this;
	}
}
