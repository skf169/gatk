package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;

import java.io.IOException;
import java.nio.file.Path;

/**
 * Tool to extract an interval from the given reference and create a new FASTA file that contains only data
 * from that interval.
 *
 * If only a contig is specified, a FASTA file containing the whole contig will be created.
 * If a contig and a start position are specified, a FASTA file containing all bases starting from the given start position to the end of the contig will be created.
 * If a contig and an end position are specified, a FASTA file containing all bases starting from the beginning of the contig to the end position will be created.
 * If a contig, start, and end position are specified, a FASTA file containing all bases from the start to the end position on the given contig will be created.
 *
 * Inputs:
 *     Contig (required)
 *     Start Position (1-based, inclusive) (optional)
 *     End Position (1-based, inclusive) (optional)
 *
 * Outputs:
 *     FASTA file
 *     FASTA index
 *     FASTA sequence dictionary
 */
public class ExtractSequence extends GATKTool {

    private static final Logger logger = LogManager.getLogger(ExtractSequence.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String CONTIG_ARG_LONG_NAME = "contig";
    public static final String START_ARG_LONG_NAME = "start";
    public static final String END_ARG_LONG_NAME = "end";

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Path to output FASTA file to which sub-sequence should be written.")
    protected Path outputFilePath;

    @Argument(
            fullName = CONTIG_ARG_LONG_NAME,
            doc = "Contig from which to extract a sub-sequence of bases."
    )
    protected String contig = null;

    @Argument(
            fullName = START_ARG_LONG_NAME,
            doc = "Start position (1-based, inclusive) from which to extract a sub-sequence of bases.",
            optional = true
    )
    protected Integer start = null;

    @Argument(
            fullName = END_ARG_LONG_NAME,
            doc = "End position (1-based, inclusive) from which to extract a sub-sequence of bases.",
            optional = true
    )
    protected Integer end = null;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void traverse() {

        // Make sure our contig exists in the reference:
        if ( reference.getSequenceDictionary().getSequence(contig) == null ) {
            throw new UserException("Error: contig does not exist in the reference file: " + contig + " : " + referenceArguments.getReferencePath());
        }

        final StringBuilder outputContigNameBuilder = new StringBuilder();
        outputContigNameBuilder.append( contig );

        // Check that we have start and end values:
        if ( start == null ) {
            // Default start to 1:
            start = 1;
        }
        else {
            // Let the user know that because of the FASTA format, we need to relabel things with new start positions:
            logger.warn("Start position will be used, but not reflected in FASTA format.  Will append interval positions to contig name in resulting file.");
            outputContigNameBuilder.append("_");
            outputContigNameBuilder.append(start);
        }

        if ( end == null ) {
            // Default end to the last position in the given contig:
            end = reference.getSequenceDictionary().getSequence(contig).getSequenceLength();

            if ( start != 1 ) {
                outputContigNameBuilder.append("-");
                outputContigNameBuilder.append(end);
            }
        }

        // Create the interval to pull out:
        final SimpleInterval interval = new SimpleInterval(contig, start, end);

        // Get the bases to write:
        final ReferenceSequence refSeq = reference.queryAndPrefetch(interval);

        // Create a FASTA writer:
        try (final FastaReferenceWriter fastaReferenceWriter = new FastaReferenceWriter(outputFilePath, true, true)) {

            // Start our contig:
            fastaReferenceWriter.startSequence(outputContigNameBuilder.toString());

            // Write our bases:
            fastaReferenceWriter.appendBases(refSeq.getBases());
        }
        catch (final IOException ex) {
            throw new GATKException("Error creating/writing to FASTA output writer: " + outputFilePath.toUri().toString(), ex);
        }
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
