package org.broadinstitute.hellbender.tools.walkers.fasta;

import com.google.common.primitives.Bytes;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Create a subset of a FASTA reference sequence
 *
 * <p>This tool creates a new reference in FASTA format consisting of only those positions or intervals
 * provided in the input data set. The output format can be partially controlled using the provided command-line
 * arguments. Specify intervals with the usual -L argument to output only the reference bases within your intervals.
 * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
 * separate fasta sequence (named numerically in order).</p>
 *
 * <h3>Input</h3>
 * <p>
 * The reference and requested intervals.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A fasta file representing the requested intervals. Each interval has a description line starting with a greater-than (">") symbol followed by sequence data.
 * The description begins with the contig name followed by the beginning position on the contig.
 * <pre>
 * For example, the fasta file for contig 1 and intervals 1:3-1:4 and 1:6-1:9
 * >1 1:3
 * AT
 * >1 1:6
 * GGGG
 * </pre>
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T FastaReferenceMaker \
 *   -R reference.fasta \
 *   -o output.fasta \
 *   -L input.intervals
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "todo",
        oneLineSummary = "todo",
        programGroup = ReferenceProgramGroup.class
)
public class FastaReferenceMaker extends ReferenceWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "path to write the fasta to")
    public String output;

    @Argument(fullName="line-width", doc="Maximum length of sequence to write per line", optional=true)
    public int basesPerLine = FastaReferenceWriter.DEFAULT_BASES_PER_LINE;

    protected FastaReferenceWriter writer;
    private int contigCount = 0;
    private int currentSequenceStartPosition = 0;
    private SimpleInterval lastPosition = null;
    private List<Byte> sequence = new ArrayList<>(10000);

    @Override
    public void onTraversalStart() {
        final Path path = IOUtils.getPath(output);
        try {
            writer = new FastaReferenceWriter(path, basesPerLine, true, true);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create " + output + ", encountered exception: " + e.getMessage(), e);
        }
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if( writer != null){
            try {
                writer.close();
            } catch (IOException | IllegalStateException e ) {
                if( e.getSuppressed().length == 0) {
                    throw new UserException("Failed to close fasta writer for " + output + " due to " + e.getMessage() + ".", e);
                }
            }
        }
    }

    @Override
    public void apply(ReferenceContext referenceContext, ReadsContext read, FeatureContext featureContext) {
        if(lastPosition == null ){
            initializeNewSequence(referenceContext);
        } else if ( !lastPosition.withinDistanceOf(referenceContext.getInterval(), 1)) {
            finalizeSequence();
            initializeNewSequence(referenceContext);
        }
        addToSequence(referenceContext);
    }

    private void addToSequence(ReferenceContext referenceContext) {
        sequence.add(referenceContext.getBase());
        lastPosition = referenceContext.getInterval();
    }

    private void finalizeSequence() {
        final String description = lastPosition.getContig() + ":" + currentSequenceStartPosition + "-" + lastPosition.getEnd();
        try {
            writer.appendSequence(String.valueOf(contigCount), description, basesPerLine, Bytes.toArray(sequence));
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Failed while writing " + output + ".", e);
        }
    }

    private void initializeNewSequence(ReferenceContext currentPosition) {
        lastPosition = currentPosition.getInterval();
        contigCount++;
        currentSequenceStartPosition = lastPosition.getStart();
        sequence.clear();
    }

    @Override
    public Object onTraversalSuccess(){
        finalizeSequence();
        return null;
    }
}