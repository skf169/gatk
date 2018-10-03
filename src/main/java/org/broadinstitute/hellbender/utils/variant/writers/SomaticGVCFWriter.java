package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.List;

/**
 * Genome-wide VCF writer for somatic (Mutect2) output
 * Merges reference blocks based on TLOD
 */
final public class SomaticGVCFWriter extends GVCFWriter {
    @VisibleForTesting
    protected int partitionPrecision;  //number of decimal places to use for TLOD block ranges

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param underlyingWriter the ultimate destination of the GVCF records
     * @param lodPartitions     a list of TLOD partitions, this list must be non-empty and every element must be larger than previous element
     */
    public SomaticGVCFWriter(final VariantContextWriter underlyingWriter, final List<Number> lodPartitions) {
        super(underlyingWriter, lodPartitions, HomoSapiensConstants.DEFAULT_PLOIDY);
    }

    public SomaticGVCFWriter(final VariantContextWriter underlyingWriter, final List<Number> lodPartitions, final int partitionPrecision) {
        super(underlyingWriter, lodPartitions, HomoSapiensConstants.DEFAULT_PLOIDY);
        this.partitionPrecision = partitionPrecision;
    }

    @VisibleForTesting
    protected int convertLODtoInt(final double LOD) {
        return (int)Math.floor(LOD * Math.pow(10, partitionPrecision));
    }

    /**
     * Create {@link HomRefBlock}s which will collectively accept variants of any genotype quality
     *
     * Each individual block covers a band of tumor LODs with the splits between bands occurring at values in {@code gqPartitions}.
     * There will be {@code gqPartitions.size() +1} bands produced
     *
     * @param gqPartitions proposed TLOD partitions as Doubles in LOD-space
     * @return a list of HomRefBlocks accepting bands of genotypes qualities split at the points specified in gqPartitions
     */
    @Override
    @VisibleForTesting
    public RangeMap<Integer,Range<Integer>> parsePartitions(final List<Number> gqPartitions) {
        calculatePartitionPrecision(gqPartitions);
        Utils.nonEmpty(gqPartitions);
        Utils.containsNoNull(gqPartitions, "The list of TLOD partitions contains a null integer");
        final RangeMap<Integer, Range<Integer>> result = TreeRangeMap.create();
        int lastThreshold = -Integer.MAX_VALUE;
        for (final Number num : gqPartitions) {
            final double value = num.doubleValue();
            final int intThreshold = convertLODtoInt(value);
            result.put(Range.closedOpen(lastThreshold, intThreshold), Range.closedOpen(lastThreshold, intThreshold));
            lastThreshold = intThreshold;
       }
        result.put(Range.closedOpen(lastThreshold, Integer.MAX_VALUE), Range.closedOpen(lastThreshold, Integer.MAX_VALUE));
        return result;
    }

    private void calculatePartitionPrecision(final List<Number> gqPartitions) {
        double smallestDelta = Double.POSITIVE_INFINITY;
        double lastLOD = Double.NEGATIVE_INFINITY;
        for (final Number num : gqPartitions) {
            final double value = num.doubleValue();
            Utils.validateArg(value != lastLOD, String.format("The value %f appears more than once in the list of TLOD partitions.", value));
            Utils.validateArg(value > lastLOD, String.format("The list of TLOD partitions is out of order. Previous value is %f but the next is %f.", lastLOD, value));
            final double delta = value - lastLOD;
            if (delta < smallestDelta) {
                smallestDelta = delta;
            }
            lastLOD = value;
        }
        partitionPrecision = (int)Math.ceil(-Math.log10(smallestDelta));
    }

    @Override
    boolean genotypeCanBeMergedInCurrentBlock(final Genotype g) {
        final TLODBlock currentTLODBlock = (TLODBlock)currentBlock;
        final double TLOD = Double.parseDouble(g.getExtendedAttribute(GATKVCFConstants.TUMOR_LOD_KEY).toString());
        return currentTLODBlock != null
                && currentTLODBlock.withinBounds(convertLODtoInt(TLOD));
    }

    /**
     * Helper function to create a new HomRefBlock from a variant context and current genotype
     *
     * @param vc the VariantContext at the site where want to start the band
     * @param g  the genotype of the sample from vc that should be used to initialize the block
     * @return a newly allocated and initialized block containing g already
     */
    @Override
    GVCFBlock createNewBlock(final VariantContext vc, final Genotype g) {
        // figure out the GQ limits to use based on the LOD of g
        final double lod = Double.parseDouble(g.getExtendedAttribute(GATKVCFConstants.TUMOR_LOD_KEY).toString());
        final Range<Integer> partition = gqPartitions.get(convertLODtoInt(lod));

        if( partition == null) {
            throw new GATKException("LOD for genotype " + g + " from " + vc + " didn't fit into any partition");
        }

        // create the block, add g to it, and return it for use
        final TLODBlock block = new TLODBlock(vc, partition.lowerEndpoint(), partition.upperEndpoint(), partitionPrecision);
        block.add(vc.getStart(), g);
        return block;
    }

}