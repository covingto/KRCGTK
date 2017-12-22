package org.bcm.hgsc.cmdline;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface CommandLineDescription {
    String value() default "No description listed";
}
