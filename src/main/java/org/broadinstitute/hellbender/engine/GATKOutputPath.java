package org.broadinstitute.hellbender.engine;

import java.io.Serializable;

/**
 * A marker class to use for GATK tool command line arguments that are outputs. These can
 * have an optional name supplied on the command line, as well as one or more optional tag/value pairs.
 */
public class GATKOutputPath extends GATKPathSpecifier implements Serializable {
    private static final long serialVersionUID = 1L;

    public GATKOutputPath(final String uriString) {
        super(uriString);
    }

}
