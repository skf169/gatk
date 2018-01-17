package org.broadinstitute.hellbender.tools.walkers;

import com.intel.genomicsdb.reader.GenomicsDBFeatureReader;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBConstants;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.getVariantContexts;
import static org.testng.Assert.*;

/**
 * Created by gauthier on 3/9/18.
 */
public class GnarlyGenotyperIntegrationTest extends CommandLineProgramTest {
    private static final VCFHeader VCF_HEADER = VariantContextTestUtils.getCompleteHeader();
    private static final String HG_00096 = largeFileTestDir + "gvcfs/HG00096.g.vcf.gz";
    private static final String HG_00268 = largeFileTestDir + "gvcfs/HG00268.g.vcf.gz";
    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();



    @DataProvider(name="VCFdata")
    public Object[][] getVCFdata() {
        return new Object[][]{
                // Simple Test, spanning deletions
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                        getTestFile("fiveSampleTest.vcf"), "chr20:250865-348163", NO_EXTRA_ARGS, b38_reference_20_21}
        };
    }


    @Test (dataProvider = "VCFdata")
    public void testUsingGenomicsDB(File[] inputs, File expected, String interval, List<String> additionalArguments, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(Arrays.asList(inputs), new SimpleInterval(interval));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        //copy in the extra files until Protobuf update in PR #4645 is ready
        Runtime.getRuntime().exec("cp " + getTestFile("vidmap.updated.json").getAbsolutePath() +" "+ tempGenomicsDB.getAbsolutePath() + "/vidmap.json");

        File output = runTool(genomicsDBUri, interval, reference, additionalArguments);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        VariantContextTestUtils.assertEqualVariants(actualVC, expectedVC);
    }

    protected File runTool(String input, String interval, String reference, List<String> additionalArguments) {
        final File output = createTempFile("GnarlyGenotyper", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
        .addArgument("V", input)
        .addArgument("L", interval);
        args.addOutput(output);

        additionalArguments.forEach(args::add);

        runCommandLine(args);
        return output;
    }

}