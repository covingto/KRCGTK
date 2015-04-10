package org.bcm.hgsc.utils;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.LogRecord;

public class Settings {
	public static int threadCount = 4;
	public static SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	public static boolean debug = false;
	
	public static Formatter defautlFormatter() {
		return new KRCGTKFormatter();
	}
	
	public static class KRCGTKFormatter extends Formatter {
		private static final DateFormat df = new SimpleDateFormat("dd/MM/yyyy hh:mm:ss.SSS");
		
		public String format(LogRecord record){
			StringBuilder builder = new StringBuilder(1000);
			builder.append(df.format(new Date(record.getMillis()))).append(": ");
			builder.append("[ ").append(record.getLevel()).append(" ] ");
			builder.append(this.formatMessage(record)).append(" ");
			builder.append(" - ").append(record.getSourceClassName()).append(".").append(record.getSourceMethodName());
			builder.append("\n");
			return builder.toString();
		}
		
		public String getHead(Handler h){
			return super.getHead(h);
		}
		
		public String getTail(Handler h){
			return super.getTail(h);
		}
	}
}
