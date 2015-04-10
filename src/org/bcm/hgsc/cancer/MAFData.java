package org.bcm.hgsc.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

public class MAFData {
	private Map<String, String> values = new HashMap<String, String>();
	
	public MAFData(List<String> h, List<String> lsplit) throws IOException {
		// TODO Auto-generated constructor stub
		
		values.put("Failed_Reason", "Pass"); // All maf records should have a failure status
		if (h.size() != lsplit.size()){
			throw new IOException("Field lengths differ!");
		}
		for (int i = 0; i < h.size(); i++){
			values.put(h.get(i), lsplit.get(i));
		}
	}
	
	public String get(String key){
		return values.get(key);
	}

	public void set(String key, String value) {
		// TODO Auto-generated method stub
		values.put(key, value);
	}
	
	public String toString(List<String> outheader){
		// return \t join of values for h in header
		List<String> tlist = new ArrayList<String>(outheader.size());
		for (String h : outheader ) {
			tlist.add(values.get(h));
		}
		return StringUtils.join(tlist, "\t");
	}

}
