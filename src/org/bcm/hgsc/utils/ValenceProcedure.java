package org.bcm.hgsc.utils;

import java.util.Map;
import java.util.Set;

import org.w3c.dom.Element;

public class ValenceProcedure extends ValenceStep{
	ValenceProcedure(Element action, String name, Map<String, String> parentReplacements, Set<String> parentDependencies){
		super(name, parentReplacements, parentDependencies);
	}
}
