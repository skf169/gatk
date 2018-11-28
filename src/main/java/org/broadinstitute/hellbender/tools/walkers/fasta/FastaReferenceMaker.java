//package org.broadinstitute.hellbender.tools.walkers.fasta;
//
//import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
//import org.broadinstitute.barclay.help.DocumentedFeature;
//import org.broadinstitute.hellbender.utils.GenomeLoc;
//import org.broadinstitute.hellbender.utils.collections.Pair;
//import org.broadinstitute.hellbender.utils.commandline.Argument;
//import org.broadinstitute.hellbender.utils.commandline.Output;
//import org.broadinstitute.hellbender.utils.contexts.AlignmentContext;
//import org.broadinstitute.hellbender.utils.contexts.ReferenceContext;
//import org.broadinstitute.hellbender.utils.refdata.RefMetaDataTracker;
//import picard.cmdline.programgroups.ReferenceProgramGroup;
//
//import java.io.PrintStream;
//
///**
// * Create a subset of a FASTA reference sequence
// *
// * <p>This tool creates a new reference in FASTA format consisting of only those positions or intervals
// * provided in the input data set. The output format can be partially controlled using the provided command-line
// * arguments. Specify intervals with the usual -L argument to output only the reference bases within your intervals.
// * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
// * separate fasta sequence (named numerically in order).</p>
// *
// * <h3>Input</h3>
// * <p>
// * The reference and requested intervals.
// * </p>
// *
// * <h3>Output</h3>
// * <p>
// * A fasta file representing the requested intervals. Each interval has a description line starting with a greater-than (">") symbol followed by sequence data.
// * The description begins with the contig name followed by the beginning position on the contig.
// * <pre>
// * For example, the fasta file for contig 1 and intervals 1:3-1:4 and 1:6-1:9
// * >1 1:3
// * AT
// * >1 1:6
// * GGGG
// * </pre>
// * </p>
// *
// * <h3>Usage example</h3>
// * <pre>
// * java -jar GenomeAnalysisTK.jar \
// *   -T FastaReferenceMaker \
// *   -R reference.fasta \
// *   -o output.fasta \
// *   -L input.intervals
// * </pre>
// *
// */
//@DocumentedFeature
//@CommandLineProgramProperties(
//  summary = "todo",
//  oneLineSummary = "todo",
//  programGroup = ReferenceProgramGroup.class
//)
//public class FastaReferenceMaker extends ReferenceWalker {
//
//  @Output PrintStream out;
//
//  @Argument(fullName="lineWidth", shortName="lw", doc="Maximum length of sequence to write per line", required=false)
//  public int fastaLineWidth=60;
//
//  /**
//   *  Please note that when using this argument adjacent intervals will automatically be merged.
//   */
//  @Argument(fullName="rawOnelineSeq", shortName="raw", doc="Print sequences with no FASTA header lines, one line per interval (i.e. lineWidth = infinity)", required=false)
//  public boolean fastaRawSeqs=false;
//
//  protected FastaSequence fasta;
//
//  public void initialize() {
//    if (fastaRawSeqs) fastaLineWidth = Integer.MAX_VALUE;
//    fasta = new FastaSequence(out, fastaLineWidth, fastaRawSeqs);
//  }
//
//  public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
//    return new Pair<GenomeLoc, String>(context.getLocation(), String.valueOf((char)ref.getBase()));
//  }
//
//  public GenomeLoc reduceInit() {
//    return null;
//  }
//
//  public GenomeLoc reduce(Pair<GenomeLoc, String> value, GenomeLoc sum) {
//    if ( value == null )
//      return sum;
//
//    // if there is no interval to the left, then this is the first one
//    if ( sum == null ) {
//      sum = value.first;
//      fasta.setName(fasta.getName() + " " + sum.toString());
//      fasta.append(value.second);
//    }
//    // if the intervals are not contiguous, print out the leftmost one and start a new one
//    // (end of contig or new interval)
//    else if ( value.first.getStart() != sum.getStop() + 1 || ! value.first.getContig().equals(sum.getContig()) ) {
//      fasta.flush();
//      sum = value.first;
//      fasta.setName(fasta.getName() + " " + sum.toString());
//      fasta.append(value.second);
//    }
//    // otherwise, merge them
//    else {
//      sum = sum.setStop(sum, value.first.getStop());
//      fasta.append(value.second);
//    }
//    return sum;
//  }
//
//  public void onTraversalDone(GenomeLoc sum) {
//    fasta.flush();
//  }
//}