package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

public class AnalyzeMITESeq extends GATKTool {
    @Argument(doc = "minimum quality score for analyzed portion of read", fullName = "min-q")
    private static int minQ = 30;

    @Argument(doc = "minimum size of analyzed portion of read", fullName = "min-length")
    private static int minLength = 15;

    @Argument(doc = "minimum number of wt calls flanking variant", fullName = "flanking-length")
    private static int flankingLength = 18;

    @Argument(doc = "reference indices of the ORF (1-based, closed), for example, '134-180,214-238'", fullName = "orf")
    private static String orfCoords;

    private byte[] refSeq;
    private List<Exon> exonList;
    private final HopscotchMap<SNVCollectionCount, Integer, SNVCollectionCount> variationCounts = new HopscotchMap<>(1000000);

    @Override
    public boolean requiresReads() { return true; }
    @Override
    public boolean requiresReference() { return true; }

    // describes an exon as a pair of offsets (0-based, half-open) on the reference sequence.
    private final static class Exon {
        private final int start;
        private final int end;

        public Exon( final int begin, final int end ) {
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
            return refIndex + ":" + (char)refCall + "->" + (char)variantCall;
        }
    }

    public static final class SNVCollectionCount
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

    @Override
    public void traverse() {
        initializeRefSeq();
        initializeExons();

        getTransformedReadStream(ReadFilterLibrary.ALLOW_ALL_READS).forEach(this::processRead);

        List<SNVCollectionCount> variationEntries = new ArrayList<>(variationCounts.size());
        variationEntries.addAll(variationCounts);
        variationEntries.sort((a,b) -> {
            int result = Integer.compare(b.getCount(), a.getCount()); // descending order of count
            if ( result == 0 ) result = a.compareTo(b);
            return result;
        });
        variationEntries.forEach(System.out::println);
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
            switch ( refSeq[idx] & 0xDF ) { // make into lower case
                case 'A': case 'C': case 'G': case 'T':
                    break;
                default:
                    throw new UserException("Reference sequence contains something other than A, C, G, and T.");
            }
        }
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
                final Exon exon = new Exon(begin-1, end);
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

        final int orfLen = exonList.stream().mapToInt(Exon::size).sum();
        if ( (orfLen % 3) != 0 ) {
            throw new UserException("ORF length must be divisible by 3.");
        }

        // it's helpful to have this 0-length sentinel at the end of the list
        exonList.add(new Exon(Integer.MAX_VALUE, Integer.MAX_VALUE));
    }

    private void processRead( final GATKRead read ) {
        // ignore unaligned reads
        if ( read.isUnmapped() ) return;

        // find pieces of the read that can be analyzed
        final byte[] quals = read.getBaseQualities();

        // advance start until we find a high quality call
        int start = 0;
        while ( start < quals.length ) {
            if ( quals[start] < minQ ) {
                start += 1;
                continue;
            }
            // advance end until we find a low quality call
            int end = start;
            while ( ++end < quals.length ) {
                if ( quals[end] < minQ ) break;
            }

            // analyze high quality pieces of sufficient length
            if ( end-start >= minLength ) {
                analyze(read, start, end);
            }

            // find another high quality region, starting where we left off
            start = end + 1;
        }
    }

    private void analyze( final GATKRead read, final int start, final int end ) {
        List<SNV> variations = new ArrayList<>();
        final byte[] readSeq = read.getBasesNoCopy();
        int readIndex = 0;
        int refIndex = read.getStart() - 1; // 0-based numbering
        final Iterator<CigarElement> cigarIterator = read.getCigar().getCigarElements().iterator();
        CigarElement cigarElement = cigarIterator.next();
        CigarOperator cigarOperator = cigarElement.getOperator();
        int cigarElementIdx = 0;
        while ( readIndex < end ) {
            if ( cigarElementIdx == cigarElement.getLength() ) {
                cigarElement = cigarIterator.next();
                cigarOperator = cigarElement.getOperator();
                cigarElementIdx = 0;
            }
            if ( readIndex >= start ) {
                if ( !cigarOperator.consumesReadBases() ) {
                    variations.add(new SNV(refIndex, refSeq[refIndex], (byte)'-'));
                } else if ( !cigarOperator.consumesReferenceBases() ) {
                    variations.add(new SNV(refIndex, (byte)'-', readSeq[readIndex]));
                } else if ( (readSeq[readIndex] & 7) != (refSeq[refIndex] & 7) ) {
                    variations.add(new SNV(refIndex, refSeq[refIndex], readSeq[readIndex]));
                }
            }
            if ( cigarOperator.consumesReadBases() ) readIndex += 1;
            if ( cigarOperator.consumesReferenceBases() ) refIndex += 1;
            cigarElementIdx += 1;
        }
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
