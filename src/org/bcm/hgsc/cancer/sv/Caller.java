package org.bcm.hgsc.cancer.sv;

import java.io.Serializable;
import java.util.UUID;

public class Caller implements Serializable{
	/**
	 * 	serial version as of 14 March 2014
	 */
	private static final long serialVersionUID = 1L;
	private final String name;
	private final int buffer;
	private final String params;
	private final String uuid;
	private final String platform;
	
	public Caller(String name, String params, int buffer, String platform){
		this.name = name;
		this.params = params;
		this.buffer = buffer;
		this.platform = platform;
		this.uuid = UUID.randomUUID().toString();
	}
	
	public Caller(String name, String params, int buffer, String platform, String uuid){
		this.name = name;
		this.params = params;
		this.buffer = buffer;
		this.uuid = uuid;
		this.platform = platform;
	}
	
	public String getPlatform(){
		return this.platform;
	}

	public String getName() {
		return name;
	}

	public int getBuffer() {
		return buffer;
	}

	public String getParams() {
		return params;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + buffer;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + ((params == null) ? 0 : params.hashCode());
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
		if (!(obj instanceof Caller)) {
			return false;
		}
		Caller other = (Caller) obj;
		if (buffer != other.buffer) {
			return false;
		}
		if (name == null) {
			if (other.name != null) {
				return false;
			}
		} else if (!name.equals(other.name)) {
			return false;
		}
		if (params == null) {
			if (other.params != null) {
				return false;
			}
		} else if (!params.equals(other.params)) {
			return false;
		}
		return true;
	}

	public String getUUID() {
		return this.uuid;
	}
}
