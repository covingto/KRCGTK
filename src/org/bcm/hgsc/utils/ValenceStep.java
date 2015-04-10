package org.bcm.hgsc.utils;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

//===========
public abstract class ValenceStep {
	protected final String name;
	protected final Set<String> dependencies = new HashSet<String>();
	protected final Map<String, String> replacements = new HashMap<String, String>();
	ValenceStep(String name, Map<String, String> parentReplacements, Set<String> parentDependencies){
		this.name = name;
		this.replacements.putAll(parentReplacements);
		this.dependencies.addAll(parentDependencies);
	}

	ValenceStep(String name){
		this.name = name;
	}

	public synchronized boolean dependenciesCleared(Set<String> cleared){
		for (String d : dependencies){
			if (!cleared.contains(d)){ return false; }
		}
		return true;
	}
}
