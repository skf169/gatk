package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "(Experimental) Processes reads from a MITESeq experiment.",
        oneLineSummary = "(EXPERIMENTAL) Processes reads from a MITESeq experiment.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public class AnalyzeMITESeq extends GATKTool {
    @Argument(doc = "minimum quality score for analyzed portion of read", fullName = "min-q")
    private static int minQ = 30;

    @Argument(doc = "minimum size of analyzed portion of read", fullName = "min-length")
    private static int minLength = 15;

    @Argument(doc = "minimum number of wt calls flanking variant", fullName = "flanking-length")
    private static int flankingLength = 18;

    @Argument(doc = "reference indices of the ORF (1-based, closed), for example, '134-180,214-238'", fullName = "orf")
    private static String orfCoords;

    @Argument(doc = "minimum number of observations of reported variants", fullName = "min-variant-obs")
    private static int minVariantObservations = 8;

    @Argument(doc = "codon translation (a string of 64 amino acid codes", fullName = "codon-translation")
    private static String codonTranslation = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF";

    @Argument(doc = "output file prefix", fullName = "output-file-prefix", shortName = "O")
    private static String outputFilePrefix;

    private byte[] refSeq;
    private List<Interval> exonList;
    private final HopscotchMap<SNVCollectionCount, Integer, SNVCollectionCount> variationCounts = new HopscotchMap<>(10000000);

    private int nReadsTotal = 0;
    private int nReadsUnmapped = 0;
    private int nReadsLowQuality = 0;
    private int nReadsWithLowQualityVariation = 0;
    private int[] refCoverage;
    private int[][] codonCounts;

    private static final int LOWERCASE_MASK = 0xDF;
    private static final int FRAME_SHIFTING_INDEL_INDEX = 64;
    private static final int FRAME_PRESERVING_INDEL_INDEX = 65;
    private static final int CODON_COUNT_ROW_SIZE = 66;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        initializeRefSeq();
        initializeExons();
    }

    @Override
    public void traverse() {
        getTransformedReadStream(ReadFilterLibrary.ALLOW_ALL_READS).forEach(this::processRead);
    }

    @Override
    public Object onTraversalSuccess() {
        long outputSize = variationCounts.stream().filter(entry -> entry.getCount() >= minVariantObservations).count();
        List<SNVCollectionCount> variationEntries = new ArrayList<>((int)outputSize);
        for ( final SNVCollectionCount entry : variationCounts ) {
            if ( entry.getCount() > minVariantObservations ) {
                variationEntries.add(entry);
            }
        }
        variationEntries.sort((a,b) -> {
            int result = Integer.compare(b.getCount(), a.getCount()); // descending order of count
            if ( result == 0 ) result = a.compareTo(b);
            return result;
        });
        final String variantsFile = outputFilePrefix + ".variantCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(variantsFile))) ) {
            for ( final SNVCollectionCount entry : variationEntries ) {
                writer.write(entry.toString());
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+variantsFile, ioe);
        }

        final String codonsFile = outputFilePrefix + ".codonCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(codonsFile))) ) {
            final int nCodons = codonCounts.length;
            for ( int idx = 0; idx != nCodons; ++idx ) {
                final int[] rowCounts = codonCounts[idx];
                writer.write(Integer.toString(idx));
                for ( final int count : rowCounts ) {
                    writer.write('\t');
                    writer.write(Integer.toString(count));
                }
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+codonsFile, ioe);
        }

        final String aaFile = outputFilePrefix + ".aaCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(codonsFile))) ) {
            final int nCodons = codonCounts.length;
            for ( int idx = 0; idx != nCodons; ++idx ) {
                final int[] rowCounts = codonCounts[idx];
                final SortedMap<Character,Integer> aaCounts = new TreeMap<>();
                for ( int codonValue = 0; codonValue != 64; ++codonValue ) {
                    aaCounts.merge(codonTranslation.charAt(codonValue), rowCounts[codonValue], Integer::sum);
                }
                writer.write(Integer.toString(idx));
                for ( final Map.Entry<Character,Integer> entry : aaCounts.entrySet() ) {
                    writer.write('\t');
                    writer.write(entry.getValue().toString());
                }
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+aaFile, ioe);
        }
        return null;
    }

    // describes an interval as a pair of offsets (0-based, half-open) on the reference sequence.
    private final static class Interval {
        private final int start;
        private final int end;

        public Interval( final int begin, final int end ) {
            this.start = begin;
            this.end = end;
        }

        public int getStart() { return start; }
        public int getEnd() { return end; }
        public int size() { return end - start; }
    }

    private static final class SNV implements Comparable<SNV> {
        private final int refIndex;
        private final byte refCall;
        private final byte variantCall;

        public SNV( final int refIndex, final byte refCall, final byte variantCall ) {
            this.refIndex = refIndex;
            this.refCall = refCall;
            this.variantCall = variantCall;
        }

        public int getRefIndex() { return refIndex; }
        public byte getRefCall() { return refCall; }
        public byte getVariantCall() { return variantCall; }

        @Override
        public int hashCode() {
            return 47*(47*(47*refIndex + refCall) + variantCall);
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNV && equals((SNV)obj);
        }

        public boolean equals( final SNV that ) {
            return this.refIndex == that.refIndex &&
                    this.refCall == that.refCall &&
                    this.variantCall == that.variantCall;
        }

        @Override
        public int compareTo( SNV that ) {
            int result = Integer.compare(this.refIndex, that.refIndex);
            if ( result == 0 ) result = Byte.compare(this.refCall, that.refCall);
            if ( result == 0 ) result = Byte.compare(this.variantCall, that.variantCall);
            return result;
        }

        @Override
        public String toString() {
            return (refIndex+1) + ":" + (char)refCall + ">" + (char)variantCall;
        }
    }

    private static final class SNVCollectionCount
            implements Map.Entry<SNVCollectionCount, Integer>, Comparable<SNVCollectionCount> {
        private static SNV[] emptyArray = new SNV[0];
        private final SNV[] snvs;
        private int count;
        private final int hash;

        public SNVCollectionCount( final List<SNV> snvs ) {
            this.snvs = snvs.toArray(emptyArray);
            this.count = 1;
            int hashVal = 0;
            for ( final SNV snv : snvs ) {
                hashVal = 47*hashVal + snv.hashCode();
            }
            hash = 47*hashVal;
        }

        @Override
        public SNVCollectionCount getKey() { return this; }

        @Override
        public Integer getValue() { return count; }

        @Override
        public Integer setValue( Integer value ) {
            Integer result = count;
            count = value;
            return result;
        }

        public int getCount() { return count; }
        public void bumpCount() { count += 1; }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNVCollectionCount && equals((SNVCollectionCount)obj);
        }

        public boolean equals( final SNVCollectionCount that ) {
            return this.hash == that.hash && Arrays.equals(this.snvs, that.snvs);
        }

        @Override
        public int hashCode() { return hash; }

        @Override
        public int compareTo( final SNVCollectionCount that ) {
            final int minSize = Math.min(this.snvs.length, that.snvs.length);
            int result = 0;
            for ( int idx = 0; idx != minSize; ++idx ) {
                result = this.snvs[idx].compareTo(that.snvs[idx]);
                if ( result != 0 ) break;
            }
            if ( result == 0 ) result = Integer.compare(this.snvs.length, that.snvs.length);
            return result;
        }

        @Override
        public String toString() {
            final StringBuilder sb = new StringBuilder(Integer.toString(count));
            String sep = "\t";
            for ( final SNV snv : snvs ) {
                sb.append(sep);
                sep = ", ";
                sb.append(snv);
            }
            return sb.toString();
        }
    }

    private static final class CodonTracker {
        private final Iterator<Interval> exonIterator;
        private Interval currentExon;
        private final int firstCodonIndex;
        private int codonPhase;
        private int codonValue;
        private int indelCount;
        private final List<Integer> codonValues;
        private int indelIndex;

        private static final int NO_INDEL_INDEX = -1;

        public CodonTracker( final List<Interval> exonList, final int refIndex ) {
            exonIterator = exonList.iterator();
            currentExon = exonIterator.next();
            int codonIndex = 0;
            while ( refIndex >= currentExon.getEnd() ) {
                codonIndex += currentExon.size();
                currentExon = exonIterator.next(); // there's a sentinel to keep us from going off the end
            }
            if ( refIndex > currentExon.getStart() ) {
                codonIndex += refIndex - currentExon.getStart();
            }
            codonPhase = codonIndex % 3;
            codonIndex /= 3;
            if ( codonPhase == 0 ) {
                firstCodonIndex = codonIndex;
            } else {
                firstCodonIndex = codonIndex + 1;
                codonPhase -= 3;
            }
            codonValue = 0;
            indelCount = 0;
            codonValues = new ArrayList<>();
            indelIndex = NO_INDEL_INDEX;
        }

        public void push( final int refIndex, final byte call ) {
            if ( refIndex == currentExon.getEnd() ) {
                currentExon = exonIterator.next();
            }
            if ( refIndex < currentExon.getStart() ) {
                return;
            }

            int callCode;
            switch ( call ) {
                case 'A': callCode = 0; break;
                case 'C': callCode = 1; break;
                case 'G': callCode = 2; break;
                case 'T': callCode = 3; break;
                case '+': indelCount += 1; return;
                case '-': indelCount -= 1; return;
                default: throw new GATKException("high quality call with value " + (char)call);
            }

            if ( indelCount != 0 ) {
                indelIndex = firstCodonIndex + codonValues.size();
                codonValues.clear(); // no other variant calling on indel-containing reads

                // kill any further calling
                while ( exonIterator.hasNext() ) {
                    currentExon = exonIterator.next();
                }
            }

            codonValue = (codonValue << 2) | callCode;

            if ( ++codonPhase == 3 ) {
                codonValues.add( codonValue & 0x3F );
                codonPhase = 0;
            }
        }

        public void report( final int[][] codonCounts ) {
            if ( indelIndex != NO_INDEL_INDEX ) {
                final int[] counts = codonCounts[indelIndex];
                final int idx = (indelCount % 3) != 0 ? FRAME_SHIFTING_INDEL_INDEX : FRAME_PRESERVING_INDEL_INDEX;
                counts[idx] += 1;
            } else {
                int idx = firstCodonIndex;
                for ( final int codonValue : codonValues ) {
                    codonCounts[idx++][codonValue] += 1;
                }
            }
        }
    }

    private void initializeRefSeq() {
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final SAMSequenceDictionary seqDict = reference.getSequenceDictionary();
        if ( seqDict.size() != 1 ) {
            throw new UserException("Expecting a reference with a single contig. " +
                    "The supplied reference has " + seqDict.size() + " contigs.");
        }
        final SAMSequenceRecord tig0 = seqDict.getSequence(0);
        final int refSeqLen = tig0.getSequenceLength();
        final SimpleInterval wholeTig = new SimpleInterval(tig0.getSequenceName(), 1, refSeqLen);
        refSeq = Arrays.copyOf(reference.queryAndPrefetch(wholeTig).getBases(),refSeqLen);
        for ( int idx = 0; idx < refSeqLen; ++idx ) {
            switch ( refSeq[idx] &= LOWERCASE_MASK ) { // make into lower case
                case 'A': case 'C': case 'G': case 'T':
                    break;
                default:
                    throw new UserException("Reference sequence contains something other than A, C, G, and T.");
            }
        }
        refCoverage = new int[refSeq.length];
    }

    private void initializeExons() {
        exonList = new ArrayList<>();
        for ( final String coordPair : orfCoords.split(",") ) {
            final String[] coords = coordPair.split("-");
            if ( coords.length != 2 ) {
                throw new UserException("Can't interpret ORF as list of pairs of coords: " + orfCoords);
            }
            try {
                final int begin = Integer.valueOf(coords[0]);
                if ( begin < 1 ) {
                    throw new UserException("Coordinates of ORF are 1-based.");
                }
                final int end = Integer.valueOf(coords[1]);
                if ( end < begin ) {
                    throw new UserException("Found ORF end coordinate less than begin: " + orfCoords);
                }
                // convert 1-based, inclusive intervals to 0-based, half-open
                final Interval exon = new Interval(begin-1, end);
                exonList.add(exon);
            }
            catch ( final NumberFormatException nfe ) {
                throw new UserException("Can't interpret ORF coords as integers: " + orfCoords);
            }
            for ( int idx = 1; idx < exonList.size(); ++idx ) {
                if ( exonList.get(idx-1).getEnd() >= exonList.get(idx).getStart() ) {
                    throw new UserException("ORF coordinates are not sorted: " + orfCoords);
                }
            }
        }

        final int orfLen = exonList.stream().mapToInt(Interval::size).sum();
        if ( (orfLen % 3) != 0 ) {
            throw new UserException("ORF length must be divisible by 3.");
        }

        codonCounts = new int[orfLen / 3][];
        for ( int idx = 0; idx != codonCounts.length; ++idx ) {
            codonCounts[idx] = new int[CODON_COUNT_ROW_SIZE];
        }

        // it's helpful to have this 0-length sentinel at the end of the list
        exonList.add(new Interval(Integer.MAX_VALUE, Integer.MAX_VALUE));
    }

    private void processRead( final GATKRead read ) {
        // ignore unaligned reads and non-primary alignments
        if ( read.isSecondaryAlignment() || read.isSupplementaryAlignment() ) return;
        nReadsTotal += 1;
        if ( read.isUnmapped() ) {
            nReadsUnmapped += 1;
            return;
        }

        final byte[] quals = read.getBaseQualitiesNoCopy();

        // find initial end-trim
        int start = 0;
        int hiQCount = 0;
        while ( start < quals.length ) {
            if ( quals[start] < minQ ) {
                hiQCount = 0;
            } else if ( ++hiQCount == minLength ) {
                break;
            }
            start += 1;
        }
        if ( start == quals.length ) {
            nReadsLowQuality += 1;
            return;
        }
        start -= minLength - 1;

        // find final end-trim
        int end = quals.length - 1;
        hiQCount = 0;
        while ( end >= 0 ) {
            if ( quals[end] < minQ ) {
                hiQCount = 0;
            } else if ( ++hiQCount == minLength ) {
                break;
            }
            end -= 1;
        }
        end += minLength;

        analyze( read, start, end );
    }

    private void analyze( final GATKRead read, final int start, final int end ) {
        final Cigar cigar = read.getCigar();

        // reads with a soft end-clip are no good unless the clip is "off the end" of the amplicon
        if ( cigar.getLastCigarElement().getOperator() == CigarOperator.S ) {
            if ( read.getEnd() != refSeq.length - 1 ) {
                nReadsWithLowQualityVariation += 1;
                return;
            }
        }
        // reads with a soft start-clip are no good unless the clip is before the beginning of the amplicon
        if ( cigar.getFirstCigarElement().getOperator() == CigarOperator.S ) {
            if ( read.getStart() > 1 ) {
                nReadsWithLowQualityVariation += 1;
                return;
            }
        }

        final Iterator<CigarElement> cigarIterator = cigar.getCigarElements().iterator();
        CigarElement cigarElement = cigarIterator.next();
        CigarOperator cigarOperator = cigarElement.getOperator();
        int cigarElementCount = cigarElement.getLength();

        final byte[] readSeq = read.getBasesNoCopy();
        final byte[] readQuals = read.getBaseQualitiesNoCopy();

        final List<SNV> variations = new ArrayList<>();

        int refIndex = read.getStart() - 1; // 0-based numbering
        int readIndex = 0;

        CodonTracker codonTracker = null;
        final List<Interval> refCoverageList = new ArrayList<>();
        Interval currentRefCoverageInterval = null;

        while ( true ) {
            if ( readIndex >= start ) {
                if ( codonTracker == null ) {
                    codonTracker = new CodonTracker(exonList, refIndex);
                }
                if ( cigarOperator == CigarOperator.D ) {
                    variations.add(new SNV(refIndex, refSeq[refIndex], (byte)'-'));
                    codonTracker.push(refIndex, (byte)'-');
                } else if ( cigarOperator == CigarOperator.I ) {
                    // low-quality variations spoil the read
                    if ( readQuals[readIndex] < minQ ) {
                        nReadsWithLowQualityVariation += 1;
                        return;
                    }
                    variations.add(new SNV(refIndex, (byte)'-', readSeq[readIndex]));
                    codonTracker.push(refIndex, (byte)'+');
                } else if ( cigarOperator == CigarOperator.M ) {
                    byte call = (byte) (readSeq[readIndex] & LOWERCASE_MASK);
                    if (call != refSeq[refIndex]) {
                        // low-quality variations spoil the read
                        if (readQuals[readIndex] < minQ) {
                            nReadsWithLowQualityVariation += 1;
                            return;
                        }
                        variations.add(new SNV(refIndex, refSeq[refIndex], readSeq[readIndex]));
                    }
                    if ( currentRefCoverageInterval == null ) {
                        currentRefCoverageInterval = new Interval(refIndex, refIndex + 1);
                    } else if ( currentRefCoverageInterval.getEnd() == refIndex ) {
                        currentRefCoverageInterval = new Interval(currentRefCoverageInterval.getStart(), refIndex + 1);
                    } else {
                        refCoverageList.add(currentRefCoverageInterval);
                        currentRefCoverageInterval = new Interval(refIndex, refIndex + 1);
                    }
                    codonTracker.push(refIndex, call);
                } else if ( cigarOperator != CigarOperator.S ) {
                    throw new GATKException("unanticipated cigar operator");
                }
            }

            if ( cigarOperator.consumesReadBases() ) {
                if ( ++readIndex == end ) {
                    break;
                }
            }

            if ( cigarOperator.consumesReferenceBases() ) {
                refIndex += 1;
            }

            if ( --cigarElementCount == 0 ) {
                cigarElement = cigarIterator.next();
                cigarOperator = cigarElement.getOperator();
                cigarElementCount = cigarElement.getLength();
            }
        }

        if ( currentRefCoverageInterval != null ) {
            refCoverageList.add(currentRefCoverageInterval);
        }

        for ( final Interval refInterval : refCoverageList ) {
            final int refIntervalEnd = refInterval.getEnd();
            for ( int idx = refInterval.getStart(); idx != refIntervalEnd; ++idx ) {
                refCoverage[idx] += 1;
            }
        }

        codonTracker.report( codonCounts );

        if ( !variations.isEmpty() &&
                refIndex - variations.get(variations.size()-1).getRefIndex() >= flankingLength &&
                variations.get(0).getRefIndex() - (read.getStart() - 1) >= flankingLength ) {
            final SNVCollectionCount newVal = new SNVCollectionCount(variations);
            final SNVCollectionCount oldVal = variationCounts.find(newVal);
            if ( oldVal != null ) oldVal.bumpCount();
            else variationCounts.add(newVal);
        }
    }
}
