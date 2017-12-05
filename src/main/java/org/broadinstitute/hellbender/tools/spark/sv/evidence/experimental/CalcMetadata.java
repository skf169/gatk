package org.broadinstitute.hellbender.tools.spark.sv.evidence.experimental;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.io.FileBackedOutputStream;
import org.apache.logging.log4j.LogManager;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.SVReadFilter;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;

/**
 * Dump some statistics about the reads.
 *
 * <p>This tool takes a file of reads as input and calculates a bag of data about them:
 * fragment length statistics by read group, mean length, coverage, partition statistics, etc. This is the first step
 * in the workflow that FindBreakpointEvidenceSpark undertakes.</p>
 * <p>This is currently a debugging tool and it is probably not generally useful to most users.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of reads aligned to reference.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A text file describing the statistics.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk CalcMetadata \
 *     -I input_reads.bam \
 *     -O statistics.txt
 * </pre>
 */
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Dump some statistics about the reads.",
        summary =
        "This tool takes a file of reads as input and calculates a bag of data about them:" +
        " fragment length statistics by read group, mean length, coverage, partition statistics, etc. This is the first step" +
        " in the workflow that FindBreakpointEvidenceSpark undertakes.",
        programGroup = StructuralVariationSparkProgramGroup.class)
public class CalcMetadata extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "output file for metadata", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "write metadata as serialized binary data, rather than as human-readable text",
            fullName = "write-as-binary", optional = true)
    private boolean writeAsBinary = true;

    @Override public boolean requiresReads() { return true; }

    @Override protected void runTool( final JavaSparkContext ctx ) {
        final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params =
                new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();
        final ReadMetadata readMetadata =
                new ReadMetadata(Collections.emptySet(), getHeaderForReads(),
                        10000, getUnfilteredReads(),
                        new SVReadFilter(params), LogManager.getLogger(CalcMetadata.class));
        if ( !writeAsBinary ) {
            ReadMetadata.writeMetadata(readMetadata, outputFile);
        } else {
            try ( final OutputStream os = new FileOutputStream(outputFile) ) {
                new ReadMetadata.Serializer().write(new Kryo(), new Output(os), readMetadata);
            } catch ( final IOException ioe ) {
                throw new UserException("Can't create output file " + outputFile, ioe);
            }
        }
    }
}