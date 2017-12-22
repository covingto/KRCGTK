package org.bcm.hgsc.cmdline;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import io.github.lukehutch.fastclasspathscanner.FastClasspathScanner;
import io.github.lukehutch.fastclasspathscanner.scanner.ScanResult;

public class KRCGTK {
    
    private static Logger log = Logger.getLogger(KRCGTK.class.getName());
        
    private static class KRCGTKCommandLineDefaults {
        
        public static final boolean COLOR_STATUS;
        
        static {
            COLOR_STATUS = getBooleanProperty("color_status", true);
        }
        
        /** Gets a string system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
        private static String getStringProperty(final String name, final String def) {
            return System.getProperty("picard.cmdline." + name, def);
        }

        /** Gets a boolean system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
        private static boolean getBooleanProperty(final String name, final boolean def) {
            final String value = getStringProperty(name, String.valueOf(def));
            return Boolean.parseBoolean(value);
        }

        /** Gets an int system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
        private static int getIntProperty(final String name, final int def) {
            final String value = getStringProperty(name, String.valueOf(def));
            return Integer.parseInt(value);
        }

        /** Gets a File system property, prefixed with "picard.cmdline." using the default if the property does not exist. */
        private static File getFileProperty(final String name, final String def) {
            final String value = getStringProperty(name, def);
            // TODO: assert that it is readable
            return (null == value) ? null : new File(value);
    }
    }
    
    private static String initializeColor(final String color) {
        if (KRCGTKCommandLineDefaults.COLOR_STATUS) return color;
        else return "";
    }
    
    /** Provides ANSI colors for the terminal output **/
    private final static String KNRM = initializeColor("\u001B[0m"); // reset
    private final static String KBLD = initializeColor("\u001B[1m"); // Bold
    private final static String KRED = initializeColor("\u001B[31m");
    private final static String KGRN = initializeColor("\u001B[32m");
    private final static String KYEL = initializeColor("\u001B[33m");
    private final static String KBLU = initializeColor("\u001B[34m");
    private final static String KMAG = initializeColor("\u001B[35m");
    private final static String KCYN = initializeColor("\u001B[36m");
    private final static String KWHT = initializeColor("\u001B[37m");
    private final static String KBLDRED = initializeColor("\u001B[1m\u001B[31m");

	public static void main(String[] args) throws IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException {
		
	    List<Class<?>> classRefs = scanCommandLineTools();
	    // System.out.println(classRefs);
	    Map<String, Class<?>> toolmap = classRefs.stream().collect(Collectors.toMap(c -> {
	        CommandLineTool a = c.getAnnotation(CommandLineTool.class);
	        return a.value();
	    }, c -> c));
	    // if there are no args or the args have -h has the first option
	    // we just print the help and exit.
	    if (args.length == 0 || args[0].equals("-h") || args[0].equals("--list-commands")) {
	        printUsage(toolmap);
	        System.exit(10);
	    }
	    
	    // what is the first argument
	    String toolname = args[0];
	    
	    if (toolmap.containsKey(toolname)) {
	        Class<?> c = toolmap.get(toolname);
	        // Create the array of Argument Types
	        Class[] argTypes = { args.getClass(), };
	        Object passedArgv[] = { Arrays.copyOfRange(args, 1, args.length) };
	        Method m = c.getMethod("main", argTypes);
	        m.invoke(null, passedArgv);
	    } else {
	        System.out.println(KBLDRED + "COMMAND NOT FOUND: " + toolname + KNRM + "\n");
	        System.out.println("Please use the following");
	        printUsage(toolmap);
	        System.exit(10);
	    }
	}
	
	private static void printUsage(Map<String, Class<?>> toolmap) {
        final StringBuilder builder = new StringBuilder();
        
        builder.append(KBLDRED + "USAGE: KRCGTK-version.jar " + " " + KGRN + "<program name>" + KBLDRED + " [-h]\n\n" + KNRM);
        builder.append(KBLDRED + "Available Programs:\n" + KNRM);
        builder.append(KWHT + "--------------------------------------------------------------------------------------\n" + KNRM);
        for (final Entry<String, Class<?>> e : toolmap.entrySet()) {
            String name = e.getKey();
            Class<?> clazz = e.getValue();
            
            CommandLineDescription d = clazz.getAnnotation(CommandLineDescription.class);
            String dd;
            if (d == null) {
                dd = "";
            } else {
                dd = d.value();
            }
            
            if (clazz.getSimpleName().length() >= 45) {
                builder.append(String.format("%s    %s    %s%s%s\n", KGRN, name, KCYN, dd, KNRM));
            } else {
                builder.append(String.format("%s    %-45s%s%s%s\n", KGRN, name, KCYN, dd, KNRM));
            }
        }
        builder.append(KWHT + "--------------------------------------------------------------------------------------\n" + KNRM);
        System.out.print(builder.toString());
        
    }

    public static List<Class<?>> scanCommandLineTools() {
        List<Class<?>> classRefs = new ArrayList<Class<?>>();
        
	    // scan for classes that have the CommandLineTool annotation
	    new FastClasspathScanner(new String[] {"org.bcm.hgsc"})
	            .matchClassesWithAnnotation(CommandLineTool.class, 
	                    // c is a class annotated with CommandLineTool
	                    c -> {
	                        // System.out.println("Found " + c.getSimpleName());
	                        classRefs.add(c);
	                    } )  
	            .scan();
	    
	    return classRefs;
	}
	
}
