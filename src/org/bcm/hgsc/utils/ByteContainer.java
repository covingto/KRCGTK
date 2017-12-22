package org.bcm.hgsc.utils;

import java.util.Arrays;

public class ByteContainer{
	public final byte[] bytes;
	
	ByteContainer(Byte[] inbytes){
		this.bytes = new byte[inbytes.length];
		for (int i = 0; i < inbytes.length; i++){
			this.bytes[i] = inbytes[i];
		}
	}
	
	ByteContainer(byte [] inbytes){
		this.bytes = inbytes.clone();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + getOuterType().hashCode();
		result = prime * result + Arrays.hashCode(bytes);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (!(obj instanceof ByteContainer)) {
			return false;
		}
		ByteContainer other = (ByteContainer) obj;
		if (!getOuterType().equals(other.getOuterType())) {
			return false;
		}
		if (!Arrays.equals(bytes, other.bytes)) {
			return false;
		}
		return true;
	}

	private ByteContainer getOuterType() {
		return this;
	}

	@Override
	public String toString() {
		return new String(this.bytes);
	}
	
	public String toString(int offset){
		byte[] newbyte = new byte[this.bytes.length];
		for (int i = 0; i < newbyte.length; i++){
			newbyte[i] = (byte) (this.bytes[i] + offset); 
		}
		return new String(newbyte);
	}
}
