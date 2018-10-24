package org.broadinstitute.hellbender.engine;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.FileSystemNotFoundException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class PathSpecifierUnitTest {

    @DataProvider
    public Object[][] getURIs() throws IOException {
        return new Object[][] {
            // input URI, expected resulting URI string, isNIO, isPath
            {"localFile.bam",                   "file://" + getCWD() + "localFile.bam", true, true},
            {"/localFile.bam",                  "file:///localFile.bam",                true, true},
            {"file:/localFile.bam",             "file:/localFile.bam",                  true, true},
            {"file:localFile.bam",              "file:localFile.bam",                   true, false}, // opaque, but not hierarchical
            {"file://localFile.bam",            "file://localFile.bam",                 true, false}, // file URLs can't have an authority ("localFile.bam")
            {"file:///localFile.bam",           "file:///localFile.bam",                true, true},  // empty authority

            {"path/to/localFile.bam",           "file://" + getCWD() + "path/to/localFile.bam", true, true},
            {"/path/to/localFile.bam",          "file:///path/to/localFile.bam",    true, true},
            {"file:path/to/localFile.bam",      "file:path/to/localFile.bam",       true, false},
            {"file:/path/to/localFile.bam",     "file:/path/to/localFile.bam",      true, true},
            {"file://path/to/localFile.bam",    "file://path/to/localFile.bam",     true, false}, // "path" looks like an authority, but won't be treated that way
            {"file:///path/to/localFile.bam",   "file:///path/to/localFile.bam",    true, true},  // empty authority

            {CloudStorageFileSystem.URI_SCHEME + ":" + "//file.bam",                CloudStorageFileSystem.URI_SCHEME + ":" + "//file.bam", true, true},
            {CloudStorageFileSystem.URI_SCHEME + ":" + "//bucket/file.bam",         CloudStorageFileSystem.URI_SCHEME + ":" + "//bucket/file.bam", true, true},
            {CloudStorageFileSystem.URI_SCHEME + ":" + "///bucket/file.bam",        CloudStorageFileSystem.URI_SCHEME + ":" + "///bucket/file.bam", true, false},
            {CloudStorageFileSystem.URI_SCHEME + ":" + "//auth/bucket/file.bam",    CloudStorageFileSystem.URI_SCHEME + ":" + "//auth/bucket/file.bam", true, true},
            {CloudStorageFileSystem.URI_SCHEME + ":" + "//hellbender/test/resources/", CloudStorageFileSystem.URI_SCHEME + ":" + "//hellbender/test/resources/", true, true},

            // java.lang.NullPointerException: Null host not permitted if default Hadoop filesystem is not HDFS.
            //{"hdfs:/file.bam",                              "hdfs:/file.bam", true, true},
            // java.lang.AssertionError: java.net.UnknownHostException: file.bam
            {"hdfs://file.bam",                             "hdfs://file.bam", true, false},
            // java.lang.NullPointerException: Null host not permitted if default Hadoop filesystem is not HDFS.
            //{"hdfs:///file.bam",                            "hdfs:///file.bam", true, true},
            // java.lang.NullPointerException: Null host not permitted if default Hadoop filesystem is not HDFS.
            //{"hdfs:/path/to/file.bam",                      "hdfs:/path/to/file.bam", true, true},
            // java.lang.AssertionError: java.net.UnknownHostException:
            {"hdfs://path/to/file.bam",                     "hdfs://path/to/file.bam", true, false},
            // java.lang.NullPointerException: Null host not permitted if default Hadoop filesystem is not HDFS.
            //{"hdfs:///path/to/file.bam",                    "hdfs:///path/to/file.bam", true, true},
            // java.lang.AssertionError: java.net.UnknownHostException: nonexistentnamenode
            {"hdfs://nonexistentnamenode/path/to/file.bam", "hdfs://nonexistentnamenode/path/to/file.bam", true, false},
            // java.lang.AssertionError: java.net.UnknownHostException: host
            {"hdfs://userinfo@host:80/path/to/file.bam",    "hdfs://userinfo@host:80/path/to/file.bam", true, false},

            // uri must have a path: jimfs:file.bam
            {"jimfs:file.bam",                              "jimfs:file.bam", true, false},
            // java.lang.AssertionError: java.net.URISyntaxException: Expected scheme-specific part at index 6: jimfs:
            {"jimfs:/file.bam",                             "jimfs:/file.bam", true, false},
            // java.lang.AssertionError: uri must have a path: jimfs://file.bam
            {"jimfs://file.bam",                            "jimfs://file.bam", true, false},
            // java.lang.AssertionError: java.net.URISyntaxException: Expected scheme-specific part at index 6: jimfs:
            {"jimfs:///file.bam",                           "jimfs:///file.bam", true, false},
            // fake_file_system_authority doesn't exist
            {"jimfs://fake_file_system_authority/file.bam", "jimfs://fake_file_system_authority/file.bam", true, false},

            // URIs with a "#" in the path part, which is valid URI syntax, but without encoding will be
            // treated as a fragment delimiter.
            {"/project/gvcf-pcr/23232_1#1/1.g.vcf.gz",      "file:///project/gvcf-pcr/23232_1%231/1.g.vcf.gz", true, true},
            {"project/gvcf-pcr/23232_1#1/1.g.vcf.gz",       "file://" + getCWD() + "project/gvcf-pcr/23232_1%231/1.g.vcf.gz", true, true},

            // URIs that are presented by the user as "file:" URLs must already be encoded (??),otherwise we'd
            // double-encode them
            {"file:project/gvcf-pcr/23232_1#1/1.g.vcf.gz",  "file:project/gvcf-pcr/23232_1#1/1.g.vcf.gz", true, false},
            {"file:/project/gvcf-pcr/23232_1#1/1.g.vcf.gz", "file:/project/gvcf-pcr/23232_1#1/1.g.vcf.gz", true, false},

            {GenomicsDBUtils.GENOMIC_DB_URI_SCHEME + "somegdb",           GenomicsDBUtils.GENOMIC_DB_URI_SCHEME + "somegdb", false, false},
            {CloudStorageFileSystem.GCS_VIEW + ":" + "//abucket/bucket",    CloudStorageFileSystem.GCS_VIEW + ":" + "//abucket/bucket", false, false},

            {"chr1:1-100", "chr1:1-100", false, false},
            {"", "file://" + getCWD(), true, true},       // an empty path is equivalent to accessing the default directory of the default file system
            {"/", "file:///", true, true},
            {"///", "file:///", true, true},

            // URI is presented already encoded/escaped and will not be altered by the URI class treatment
            {"file:///project/gvcf-pcr/23232_1%231/1.g.vcf.g", "file:///project/gvcf-pcr/23232_1%231/1.g.vcf.g", true, true}
        };
    }

    @Test(dataProvider = "getURIs")
    public void testURIValid(final String uriString, final String expectedURIString, final boolean isNIO, final boolean isPath) {
        final PathURI pathURI = new PathSpecifier(uriString);
        Assert.assertNotNull(pathURI);
        Assert.assertEquals(pathURI.getURI().toString(), expectedURIString);
    }

    @DataProvider
    public Object[][] getInvalidURIs() {
        return new Object[][] {
                {"file://^"},
                {"file://"},
                {"underbar_is_invalid_in_scheme:///foobar"},
        };
    }

    @Test(dataProvider = "getInvalidURIs", expectedExceptions = IllegalArgumentException.class)
    public void testURIInvalid(final String invalidURIString) {
        new PathSpecifier(invalidURIString);
    }

    @Test(dataProvider = "getURIs")
    public void testIsNIOValid(final String uriString, final String expectedURIString, final boolean isNIO, final boolean isPath) {
        final PathURI pathURI = new PathSpecifier(uriString);
        Assert.assertEquals(pathURI.isNIO(), isNIO);
    }

    @DataProvider
    public Object[][] getInvalidNIOURIs() {
        return new Object[][]{
                // URIs with schemes that don't have an NIO provider
                {"unknownscheme://foobar"},
                {CloudStorageFileSystem.GCS_VIEW + ":" + "//abucket/bucket"},
                {GenomicsDBUtils.GENOMIC_DB_URI_SCHEME + "adb"},
        };
    }

    @Test(dataProvider = "getInvalidNIOURIs")
    public void testIsNIOInvalid(final String invalidNIOURI) {
        final PathURI pathURI = new PathSpecifier(invalidNIOURI);
        Assert.assertEquals(pathURI.isNIO(), false);
    }

    @Test(dataProvider = "getURIs")
    public void testIsPathValid(final String uriString, final String expectedURIString, final boolean isNIO, final boolean isPath) {
        final PathURI pathURI = new PathSpecifier(uriString);
        if (isPath) {
            Assert.assertEquals(pathURI.isPath(), isPath, pathURI.getToPathFailureReason());
        } else {
            Assert.assertEquals(pathURI.isPath(), isPath);
        }
    }

    @Test(dataProvider = "getURIs")
    public void testToPathValid(final String uriString, final String expectedURIString, final boolean isNIO, final boolean isPath) {
        final PathURI pathURI = new PathSpecifier(uriString);
        if (isPath) {
            final Path path = pathURI.toPath();
            Assert.assertEquals(path != null, isPath, pathURI.getToPathFailureReason());
        } else {
            Assert.assertEquals(pathURI.isPath(), isPath);
        }
    }

    @DataProvider
    public Object[][] getInvalidPaths() {
        return new Object[][]{
                // valid URIs that are not valid as a path

                {"file:/project/gvcf-pcr/23232_1#1/1.g.vcf.gz"},    // not encoded
                {"file://path/to/file.bam"},  // Paths.get throws IllegalArgumentException, 2 leading slashes (vs. 3)
                                              // causes "path" to be interpreted as an invalid authority name
                {"file:project/gvcf-pcr/23232_1#1/1.g.vcf.gz"},     // scheme-specific part is not hierarchichal

                // The hadoop file system provider explicitly throws an NPE if no host is specified and HDFS is not
                // the default file system
                //{"hdfs://nonexistent_authority/path/to/file.bam"},  // unknown authority "nonexistent_authority"
                {"hdfs://userinfo@host:80/path/to/file.bam"},           // UnknownHostException "host"

                // TODO NOTE: the URIs from here down are accepted by IOUtils (public static Path getPath(String uriString))
                // as valid Paths, even though they are unresolvable and otherwise pretty much useless.
                {"unknownscheme://foobar"},
                {GenomicsDBUtils.GENOMIC_DB_URI_SCHEME + "adb"},
                {CloudStorageFileSystem.GCS_VIEW + ":" + "//abucket/bucket"},

                // URIs with schemes that are backed by an valid NIO provider, but for which the
                // scheme-specific part is not valid.
                {"file://nonexistent_authority/path/to/file.bam"},  // unknown authority "nonexistent_authority"
                {"file://path/to/file.bam"},                        // unknown authority "path"
        };
    }

    @Test(dataProvider = "getInvalidPaths")
    public void testIsPathInvalid(final String invalidPathString) {
        final PathURI pathURI = new PathSpecifier(invalidPathString);
        Assert.assertFalse(pathURI.isPath());
    }

    @Test(dataProvider = "getInvalidPaths", expectedExceptions = {IllegalArgumentException.class, FileSystemNotFoundException.class})
    public void testToPathInvalid(final String invalidPathString) {
        final PathURI pathURI = new PathSpecifier(invalidPathString);
        pathURI.toPath();
    }

    @DataProvider
    public Object[][] getInputStreamURIs() throws IOException {
        return new Object[][]{
                // URIs that can be resolved to an actual test file
                {"src/test/resources/org/broadinstitute/hellbender/tools/testTextFile.txt", "Test file."},
                {"file:///" + getCWD() + "src/test/resources/org/broadinstitute/hellbender/tools/testTextFile.txt", "Test file."},

                // URIs with an embedded fragemnt delimiter ("#"); if the file scheme is included, the rest of the path
                // must already be escaped; ifno file scheme is included, the path will be escaped by the URI class
                {"src/test/resources/org/broadinstitute/hellbender/tools/testDirWith#InName/testTextFile.txt", "Test file."},
                {"file:///" + getCWD() + "src/test/resources/org/broadinstitute/hellbender/tools/testDirWith%23InName/testTextFile.txt", "Test file."},
        };
    }

    @Test(dataProvider = "getInputStreamURIs")
    public void testGetInputStream(final String uriString, final String expectedFileContents) throws IOException {
        final PathURI pathURI = new PathSpecifier(uriString);

        try (final InputStream is = pathURI.getInputStream();
             final DataInputStream dis = new DataInputStream(is)) {
            final byte[] actualFileContents = new byte[expectedFileContents.length()];
            dis.readFully(actualFileContents);

            Assert.assertEquals(new String(actualFileContents), expectedFileContents);
        }
    }

    @DataProvider
    public Object[][] getOutputStreamURIs() throws IOException {
        return new Object[][]{
                // output URIs that can be resolved to an actual test file
                { IOUtils.createTempFile("testOutputStream", ".txt").getPath()},
                {"file://" + IOUtils.createTempFile("testOutputStream", ".txt").getAbsolutePath()},
        };
    }

    @Test(dataProvider = "getOutputStreamURIs")
    public void testGetOutputStream(final String tempURIString) throws IOException {
        testStreamRoundTrip(tempURIString);
    }

    @Test
    public void testAlternateFileSystem() throws IOException {
        // create a jimfs file system and round trip through PathSpecifier/stream
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path outputPath = jimfs.getPath("alternateFileSystemTest.txt");
            testStreamRoundTrip(outputPath.toUri().toString());
        }
    }

    // Temporary test for comparison of current IOUtils implementation with new implementation
    // Disabled because IOUtils accepts several schemes as valid paths
    @Test(dataProvider = "getURIs")
    public void testCompareValidPathIOUtils(final String uriString, final String expectedURIString, final boolean isNIO, final boolean isPath) {
        final PathURI pathURI = new PathSpecifier(uriString);
        if (pathURI.isPath()) {
            Assert.assertNotNull(IOUtils.getPath(uriString));
        } else {
            final String scheme = pathURI.getScheme();
            if (scheme.equals("gendb") || scheme.equals("chr1") || scheme.equals("gcs")) {
                // IOUtils creates paths for these, but it probably shouldn't
                IOUtils.getPath(uriString);
            } else {
                //these are not "Path"-able, so don't even try IOUtils
            }
        }
    }

    private void testStreamRoundTrip(final String tempURIString) throws IOException {        final String expectedFileContents = "Test contents";
        final PathURI pathURI = new PathSpecifier(tempURIString);
        try (final OutputStream os = pathURI.getOutputStream();
             final DataOutputStream dos = new DataOutputStream(os)) {
            dos.write(expectedFileContents.getBytes());
        }

        // read it back in and make sure it matches expected contents
        try (final  InputStream is = pathURI.getInputStream();
             final DataInputStream dis = new DataInputStream(is)) {
            final byte[] actualFileContents = new byte[expectedFileContents.length()];
            dis.readFully(actualFileContents);

            Assert.assertEquals(new String(actualFileContents), expectedFileContents);
        }
    }

    private String getCWD() throws IOException {
        final File cwd = new File(".");
        return cwd.getCanonicalPath()  +"/";
    }
}
