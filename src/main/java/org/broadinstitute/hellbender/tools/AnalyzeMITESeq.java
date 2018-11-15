package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class AnalyzeMITESeq extends GATKTool {
    @Argument(doc = "minimum quality score for analyzed portion of read",
            fullName = "min-q")
    private static int minQ = 30;

    @Argument(doc = "minimum size of analyzed portion of read",
            fullName = "min-length")
    private static int minLength = 15;

    @Override
    public boolean requiresReads() { return true; }
    @Override
    public boolean requiresReference() { return true; }

    @Override
    public void traverse() {
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final SAMSequenceDictionary seqDict = reference.getSequenceDictionary();
        if ( seqDict.size() != 1 ) {
            throw new UserException("Expecting a reference with a single contig. " +
                                    "The supplied reference has " + seqDict.size());
        }
        final SAMSequenceRecord tig0 = seqDict.getSequence(0);
        final SimpleInterval wholeTig = new SimpleInterval(tig0.getSequenceName(), 1, tig0.getSequenceLength());
        final byte[] refSeq = reference.queryAndPrefetch(wholeTig).getBases();

        getTransformedReadStream(ReadFilterLibrary.ALLOW_ALL_READS).forEach(read -> processRead(read, refSeq));
    }

    private void processRead( final GATKRead read, final byte[] refSeq ) {
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
                analyze(read, start, end, refSeq);
            }

            // find another high quality region, starting where we left off
            start = end + 1;
        }
    }

    private void analyze( final GATKRead read, final int start, final int end, byte[] refSeq ) {
        for ( final CigarElement element : read.getCigar().getCigarElements() ) {
        }
    }

}
