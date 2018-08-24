package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.*;

@Test(groups = "Spark")
public class RevertSamSparkIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getToolTestDataDir() {
        return "src/test/resources/org/broadinstitute/hellbender/tools/spark/RevertSamSpark";
    }

    public static List<String> defaultAttributesToClearPlusXT = new ArrayList<String>() {{
        add(SAMTag.NM.name());
        add(SAMTag.UQ.name());
        add(SAMTag.PG.name());
        add(SAMTag.MD.name());
        add(SAMTag.MQ.name());
        add(SAMTag.SA.name()); // Supplementary alignment metadata
        add(SAMTag.MC.name()); // Mate Cigar
        add(SAMTag.AS.name());
        add("XT");
    }};

    private final File basicSamToRevert = getTestFile("revert_sam_basic.sam");
    private final File sampleLibraryOverrideSam = getTestFile("revert_sam_sample_library_override.sam");
    private final File validOutputMap = getTestFile("revert_sam_valid_output_map.txt");
    private final File nonExistentOutputMap = getTestFile("revert_sam_does_not_exist.txt");
    private final File badHeaderOutputMap = getTestFile("revert_sam_bad_header_output_map.txt");
    private final File referenceFasta = getTestFile("test.fasta");
    private final File singleEndSamToRevert = getTestFile("revert_sam_single_end.sam");
    private final File missingRGInfo = getTestFile("missing-rg-info.sam");

    private static final String revertedQualities  =
            "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";

    private static final String unmappedRead = "both_reads_present_only_first_aligns/2";



    @DataProvider(name="positiveTestData")
    public Object[][] getPostitiveTestData() {
        return new Object[][] {
                {null, false, false, true, true, null, null, Collections.EMPTY_LIST},
                {SAMFileHeader.SortOrder.queryname, false, false, true, false, "Hey,Dad!", null, defaultAttributesToClearPlusXT},
                {null, true, false, false, false, "Hey,Dad!", "NewLibraryName", defaultAttributesToClearPlusXT},
                {null, true, true, false, false, null, null, Collections.EMPTY_LIST}
        };
    }

    @Test(dataProvider="positiveTestData")
    public void basicPositiveTests(final SAMFileHeader.SortOrder so, final boolean removeDuplicates, final boolean removeAlignmentInfo,
                                   final boolean restoreOriginalQualities, final boolean outputByReadGroup, final String sample, final String library,
                                   final List<String> attributesToClear) throws Exception {

        final File output = outputByReadGroup?Files.createTempDirectory("picardRevertSamSparkTest").toFile():File.createTempFile("reverted", ".sam");
        File output0 = new File(output.getPath()+"/0.sam");
        File output1 = new File(output.getPath()+"/1.sam");
        File output2 = new File(output.getPath()+"/2.sam");
        File output3 = new File(output.getPath()+"/3.sam");
        output.deleteOnExit();
        final RevertSamSpark reverter = new RevertSamSpark();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(basicSamToRevert);
        args.addOutput(output);

        if (outputByReadGroup) {
            args.addPositionalArgument("--"+RevertSamSpark.OUTPUT_BY_READGROUP_ARG);
        }
        if (so != null) {
            args.addArgument("sort-order",so.name()); //TODO decide on sort order outputing
        }
        if (!removeAlignmentInfo) {
            args.addPositionalArgument("--"+RevertSamSpark.DONT_REMOVE_ALIGNMENT_INFORMATION_ARG);
        }
        if (sample != null) {
            args.addArgument("sample-alias",sample);
        }
        if (library != null) {
            args.addArgument("library-name",library);
        }
        for (final String attr : attributesToClear) {
            args.addArgument("attributes-to-clear",attr);
        }

        runCommandLine(args);

        if (outputByReadGroup) {
            verifyPositiveResults(output0, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "0", 2, sample, library);
            verifyPositiveResults(output1, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "1", 4, sample, library);
            verifyPositiveResults(output2, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "2", 2, sample, library);
            verifyPositiveResults(output3, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "3", 0, sample, library);
        } else {
            verifyPositiveResults(output, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, null, 8, sample, library);
        }
    }

    @Test
    public void testOutputByReadGroupWithOutputMap() throws Exception {
        final File outputDir = createTempDir("testOutputByReadGroupWithOutputMap");
        outputDir.deleteOnExit();
        // Create the output map
        final File outputMapFile = Files.createTempFile("picardRevertSamSparkTestOutputMap", ".txt").toFile();
        final PrintWriter mapWriter = new PrintWriter(outputMapFile);
        final String outputPath0 = outputDir + "/my_rg0.sam";
        final String outputPath1 = outputDir + "/rg1.sam";
        final String outputPath2 = outputDir + "/my_rg2.bam";
        final String outputPath3 = outputDir + "/my_rg3.sam";//TODO not used?
        mapWriter.println("READ_GROUP_ID\tOUTPUT");
        mapWriter.println("0\t" + outputPath0);
        mapWriter.println("2\t" + outputPath2);
        mapWriter.println("1\t" + outputPath1);
        mapWriter.println("3\t" + outputPath3);
        System.out.println("outputFile: " + outputPath0);
        System.out.println("outputFile: " + outputPath1);
        System.out.println("outputFile: " + outputPath2);
        System.out.println("outputFile: " + outputPath3);
        mapWriter.close();
        outputMapFile.deleteOnExit();

        final RevertSamSpark reverter = new RevertSamSpark();

        final String args[] = new String[] {
                "-I",basicSamToRevert.getPath(),
                "--output-by-readgroup",
                "--output-map",outputMapFile.getPath(),
                "-R",referenceFasta.getPath(),
                "--sort-order",SAMFileHeader.SortOrder.queryname.name(),
                "--"+RevertSamSpark.SAMPLE_ALIAS_ARG,"test_sample_1",
                "--"+RevertSamSpark.LIBRARY_NAME_ARG,"test_library_1",
                "--"+RevertSamSpark.ATTRIBUTE_TO_CLEAR_ARG,SAMTag.NM.name()
        };

        runCommandLine(args);

        final File output0 = new File(outputPath0);
        final File output1 = new File(outputPath1);
        final File output2 = new File(outputPath2);
        verifyPositiveResults(output0, reverter, true, true, true, true, "0", 2, "test_sample_1", "test_library_1");
        verifyPositiveResults(output1, reverter, true, true, true, true, "1", 4, "test_sample_1", "test_library_1");
        verifyPositiveResults(output2, reverter, true, true, true, true, "2", 2, "test_sample_1", "test_library_1");
    }

    @Test (expectedExceptions = UserException.class)
    public void testSingleEndSanitize() throws Exception {
        final File output = File.createTempFile("single_end_reverted", ".sam");
        output.deleteOnExit();
        final String args[] = { "-I " + singleEndSamToRevert, "-O " + output.getAbsolutePath(), "--sanitize"};
        runCommandLine(args);
    }

    private void verifyPositiveResults(
            final File outputFile,
            final RevertSamSpark reverter,
            final boolean removeDuplicates,
            final boolean removeAlignmentInfo,
            final boolean restoreOriginalQualities,
            final boolean outputByReadGroup,
            final String readGroupId,
            final int numReadsExpected,
            final String sample,
            final String library) {

        outputFile.deleteOnExit();
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(outputFile);
        final SAMFileHeader header = reader.getFileHeader();
        Assert.assertEquals(header.getSortOrder(), SAMFileHeader.SortOrder.queryname);
        Assert.assertEquals(header.getProgramRecords().size(), removeAlignmentInfo ? 0 : 1);
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (outputByReadGroup) {
            Assert.assertEquals(readGroups.size(), 1);
            Assert.assertEquals(readGroups.get(0).getId(), readGroupId);
        }
        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            if (sample != null) {
                Assert.assertEquals(rg.getSample(), sample);
            } else {
                Assert.assertEquals(rg.getSample(), "Hi,Mom!");
            }
            if (library != null) {
                Assert.assertEquals(rg.getLibrary(), library);
            } else {
                Assert.assertEquals(rg.getLibrary(), "my-library");
            }
        }
        int numReads = 0;
        for (final SAMRecord rec : reader) {
            numReads++;
            if (removeDuplicates) {
                Assert.assertFalse(rec.getDuplicateReadFlag(),
                        "Duplicates should have been removed: " + rec.getReadName());
            }

            if (removeAlignmentInfo) {
                Assert.assertTrue(rec.getReadUnmappedFlag(),
                        "Alignment info should have been removed: " + rec.getReadName());
            }

            if (restoreOriginalQualities && !unmappedRead.equals(
                    rec.getReadName() + "/" + (rec.getFirstOfPairFlag() ? "1" : "2"))) {

                Assert.assertEquals(rec.getBaseQualityString(), revertedQualities);
            } else {
                Assert.assertNotSame(rec.getBaseQualityString(), revertedQualities);
            }

            for (final SAMRecord.SAMTagAndValue attr : rec.getAttributes()) {
                if (removeAlignmentInfo || (!attr.tag.equals("PG") && !attr.tag.equals("NM")
                        && !attr.tag.equals(SAMTag.MQ.toString()))) {
                    Assert.assertFalse(reverter.attributesToClear.contains(attr.tag),
                            attr.tag + " should have been cleared.");
                }
            }
        }
        Assert.assertEquals(numReads, numReadsExpected);
        CloserUtil.close(reader);
    }

    @Test
    public void testSanitizeAndDeduplicateRecords() throws Exception {
        final File input  = File.createTempFile("test-input-santize-and-deduplicate-records", ".sam");
        final File output = File.createTempFile("test-output-santize-and-deduplicate-records", ".sam");

        // Create a SAM file that has duplicate records
        final SamReader reader = SamReaderFactory.makeDefault().open(basicSamToRevert);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, input);
        int numDuplicated = 0;
        for (final SAMRecord rec : reader) {
            writer.addAlignment(rec);
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                writer.addAlignment(rec);
                numDuplicated++;
            }
        }
        reader.close();
        writer.close();

        // Make sure some records are duplicated
        Assert.assertTrue(numDuplicated > 0);

        final String [] args = new String[]{
                "--input", input.getAbsolutePath(),
                "--sanitize",
                "--keep-first-duplicate",
                "--"+RevertSamSpark.DONT_RESTORE_ORIGINAL_QUALITIES_ARG,
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);
        verifyPositiveResults(output, new RevertSamSpark(), false, true, false, false, null, 8, null, null);
    }

    @Test(dataProvider="overrideTestData", expectedExceptions = {GATKException.class})
    public void testSampleLibraryOverride(final String sample, final String library) throws Exception {
        final File output = createTempFile("bad", ".sam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(sampleLibraryOverrideSam);
        args.addOutput(output);
        if (sample != null) {
            args.addArgument(RevertSamSpark.SAMPLE_ALIAS_ARG,sample);
        }
        if (library != null) {
            args.addArgument(RevertSamSpark.LIBRARY_NAME_ARG,library);
        }
        runCommandLine(args);
    }

    @DataProvider(name="overrideTestData")
    public Object[][] getNegativeTestData() {
        return new Object[][] {
                {"NewSample", null},
                {null, "NewLibrary"},
                {"NewSample", "NewLibrary"}
        };
    }

    @Test
    public void testValidateOutputParamsByReadGroupMapValid() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsByReadGroup(null, validOutputMap, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsByReadGroupMissingMap() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsByReadGroup(null, nonExistentOutputMap, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Cannot read"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupBadHeaderMap() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsByReadGroup(null, badHeaderOutputMap, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Invalid header"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupNoMapOrDir() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsByReadGroup(null, null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("Must provide either"), true);
    }

    @Test
    public void testValidateOutputParamsByReadGroupDirValid() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsByReadGroup(createTempDir("testValidateOutputParamsNotByReadGroupValid"), null, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupValid() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsNotByReadGroup(createTempFile("testValidateOutputParamsNotByReadGroupValid",""), null, errors);
        Assert.assertEquals(errors.size(), 0);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupNoOutput() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsNotByReadGroup(null, null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("output is required"), true);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupMap() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsNotByReadGroup(null, validOutputMap, errors);
        Assert.assertEquals(errors.size(), 2);
        Assert.assertEquals(errors.get(0).contains("Cannot provide outputMap"), true);
        Assert.assertEquals(errors.get(1).contains("output is required"), true);
    }

    @Test
    public void testValidateOutputParamsNotByReadGroupDir() {
        final List<String> errors = new ArrayList<String>();
        RevertSamSpark.ValidationUtil.validateOutputParamsNotByReadGroup(createTempDir("testValidateOutputParamsNotByReadGroupDir"), null, errors);
        Assert.assertEquals(errors.size(), 1);
        Assert.assertEquals(errors.get(0).contains("should not be a directory"), true);
    }

    @Test
    public void testAssertAllReadGroupsMappedSuccess() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");

        final Map<String, File> outputMap = new HashMap<String, File>();
        outputMap.put("rg1", new File("/foo/bar/rg1.bam"));
        outputMap.put("rg2", new File("/foo/bar/rg2.bam"));
        RevertSamSpark.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1, rg2));
        RevertSamSpark.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1));
        RevertSamSpark.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg2));
    }

    @Test(expectedExceptions = {GATKException.class})
    public void testAssertAllReadGroupsMappedFailure() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        final SAMReadGroupRecord rg3 = new SAMReadGroupRecord("rg3");

        final Map<String, File> outputMap = new HashMap<String, File>();
        outputMap.put("rg1", new File("/foo/bar/rg1.bam"));
        outputMap.put("rg2", new File("/foo/bar/rg2.bam"));
        RevertSamSpark.ValidationUtil.assertAllReadGroupsMapped(outputMap, Arrays.asList(rg1, rg2, rg3));
    }

    @Test
    public void testIsOutputMapHeaderValid() {
        boolean isValid = RevertSamSpark.ValidationUtil.isOutputMapHeaderValid(Arrays.asList("READ_GROUP_ID", "OUTPUT"));
        Assert.assertEquals(isValid, true);

        isValid = RevertSamSpark.ValidationUtil.isOutputMapHeaderValid(Arrays.asList("OUTPUT"));
        Assert.assertEquals(isValid, false);

        isValid = RevertSamSpark.ValidationUtil.isOutputMapHeaderValid(Collections.EMPTY_LIST);
        Assert.assertEquals(isValid, false);
    }

    @Test
    public void testFilePathsWithoutMapFile() {
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");

        final Map<String, File> outputMap = RevertSamSpark.getOutputMap(null, new File("/foo/bar"), ".bam", Arrays.asList(rg1, rg2), true);
        Assert.assertEquals(outputMap.get("rg1"), new File("/foo/bar/rg1.bam"));
        Assert.assertEquals(outputMap.get("rg2"), new File("/foo/bar/rg2.bam"));
    }

    @Test
    public void testFilePathsWithMapFile() {
        final Map<String, File> outputMap = RevertSamSpark.getOutputMap(validOutputMap, null, ".bam", Collections.emptyList(), true);
        Assert.assertEquals(outputMap.get("rg1"), new File("/path/to/my_rg_1.ubam"));
        Assert.assertEquals(outputMap.get("rg2"), new File("/path/to/my_rg_2.ubam"));
    }

    @Test
    public void testGetDefaultExtension() {
        Assert.assertEquals(RevertSamSpark.getDefaultExtension("this.is.a.sam", RevertSamSpark.FileType.dynamic), ".sam");
        //Assert.assertEquals(RevertSamSpark.getDefaultExtension("this.is.a.cram", RevertSamSpark.FileType.dynamic), ".cram");
        Assert.assertEquals(RevertSamSpark.getDefaultExtension("this.is.a.bam", RevertSamSpark.FileType.dynamic), ".bam");
        Assert.assertEquals(RevertSamSpark.getDefaultExtension("foo", RevertSamSpark.FileType.dynamic), ".bam");
    }

    @Test
    public void testNoRgInfoSanitize() throws Exception {
        final File output = File.createTempFile("no-rg-reverted", ".sam");
        output.deleteOnExit();
        final String [] args = new String[]{
                "-I",missingRGInfo.getAbsolutePath(),
                "--sanitize",
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);
        verifyPositiveResults(output, new RevertSamSpark(), true,  true, false, false, null, 240, null, null);
    }

}