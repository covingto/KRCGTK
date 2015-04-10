package org.bcm.hgsc.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import org.w3c.dom.DOMException;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

//===========
	class ValenceAction extends ValenceStep implements Runnable{
		private final File command;
		ValenceAction(Element action, String name, Map<String, String> parentReplacements, Set<String> parentDependencies) throws DOMException, IOException{
			super(name, parentReplacements, parentDependencies);
			// handle the dependencies
			NodeList del = action.getElementsByTagName("dependencies");
			for (int i = 0; i < del.getLength(); i++){
				Node node = del.item(i);
				if (node.getNodeType() == Node.ELEMENT_NODE){
					Element dependenciesElement = (Element) node;
					NodeList dl = dependenciesElement.getElementsByTagName("depend");
					for (int d = 0; d < dl.getLength(); d++){
						Node dn = dl.item(d);
						dependencies.add(dn.getTextContent());
					}
				}
			}
			// collect the replacements
			Element replace = ValenceExecutor.getFirstElement(action, "replacements");
			NodeList rel = replace.getChildNodes();
			for (int i = 0; i < rel.getLength(); i++){
				Node node = rel.item(i);
				if (node.getNodeType() == Node.ELEMENT_NODE){
					// all of these will be replacements
					replacements.put("${" + node.getNodeName() + "}", node.getTextContent());
				}
			}
			// handle the command
			//		This involves first writing the commands out to a file, this really isn't the best way to do this,
			//		but seems to be required since we can't actually control pipes, and redirects here
			this.command = new File(ValenceExecutor.scriptDir, name + ".sh");
			FileWriter fw = new FileWriter(this.command);
			BufferedWriter bw = new BufferedWriter(fw);
			NodeList cmdl = action.getElementsByTagName("command");
			for (int i = 0; i < cmdl.getLength(); i++){
				Node node = cmdl.item(i);
				if (node.getNodeType() == Node.ELEMENT_NODE){
					bw.write(handleReplacements(node.getTextContent()));
					bw.newLine();
				}
			}
			bw.flush();
			bw.close();
		}

		private String handleReplacements(String s){
			for (String k : replacements.keySet()){
				s = s.replaceAll(k, replacements.get(k));
			}
			return s;
		}

		@Override
		public void run() {
			Runtime rt = Runtime.getRuntime();
			System.out.println("Execing; sh " + this.command);
			try {
				Process proc = rt.exec("sh " + this.command);
				// any error message?
				StreamGobbler errorGobbler = new 
						StreamGobbler(proc.getErrorStream(), "ERROR");            

				// any output?
				StreamGobbler outputGobbler = new 
						StreamGobbler(proc.getInputStream(), "OUTPUT");

				// kick them off
				errorGobbler.start();
				outputGobbler.start();

				// any error???
				int exitVal = proc.waitFor();
				ValenceExecutor.writer.writeStreamGobbler(outputGobbler, errorGobbler, this.name + " exit value; " + exitVal);
				ValenceExecutor.complete.add(name);
				ValenceExecutor.poll();
			} catch (IOException ioe) {
				ioe.printStackTrace(); 
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

