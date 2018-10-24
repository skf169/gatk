package org.broadinstitute.hellbender.engine;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import org.broadinstitute.barclay.argparser.TaggedArgument;
import org.broadinstitute.barclay.argparser.TaggedArgumentParser;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.Objects;

/**
 * A marker class to use for GATK tool command line arguments that are input or output resources. These can
 * have an optional name supplied on the command line, as well as one or more optional tag/value pairs.
 *
 */
// TODO: need a way to tag GCS-enabled args
// TODO: should we do anything with query params (i.e. -> tags) ?
public class GATKPathSpecifier extends PathSpecifier implements TaggedArgument, Serializable {
    private static final long serialVersionUID = 1L;

    public static final String GCS_SCHEME = "gs";
    public static final String HDFS_SCHEME = "hdfs";

    private String tagName;
    private Map<String, String> tagAttributes;

    public GATKPathSpecifier(final String uriString) {
        super(uriString);
    }

    @Override
    public Path toPath() {
        // special case GCS, in case the filesystem provider wasn't installed properly but is available.
        if (CloudStorageFileSystem.URI_SCHEME.equals(getURI().getScheme())) {
            final Path tempPath = BucketUtils.getPathOnGcs(getURIString());
            setCachedPath(tempPath);
            return tempPath;
        } else {
            return super.toPath();
        }
    }

    //TODO: should this wrap the stream in a .gz in a BlockCompressedInputStream like IOUtils does?
    @Override
    public InputStream getInputStream() {
        if (!isPath()) {
            throw new UserException(getToPathFailureReason());
        }

        try {
            InputStream inputStream;
            if (getURI().getScheme().equals(GCS_SCHEME)) {
                Path p = BucketUtils.getPathOnGcs(getURIString());
                inputStream = Files.newInputStream(p);
            } else if (getURI().getScheme().equals(HDFS_SCHEME)) {
                org.apache.hadoop.fs.Path file = new org.apache.hadoop.fs.Path(getURIString());
                org.apache.hadoop.fs.FileSystem fs = file.getFileSystem(new org.apache.hadoop.conf.Configuration());
                inputStream = fs.open(file);
            } else {
                inputStream = super.getInputStream();
            }

            //TODO: should this wrap the stream in a .gz in a BlockCompressedInputStream ?
            return inputStream;
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(getRawInputString(), e);
        }
    }

    /**
     * Open a path for writing regardless of whether it's on GCS, HDFS or local disk.
     * For writing to GCS it'll use the application/octet-stream MIME type.
     *
     * @return an OutputStream that writes to the path specified by this UR.
     */
    @Override
    public OutputStream getOutputStream() {
        if (!isPath()) {
            throw new UserException(getToPathFailureReason());
        }

        try {
            OutputStream outputStream;
            if (getURI().getScheme().equals(GCS_SCHEME)) {
                Path p = BucketUtils.getPathOnGcs(getURIString());
                outputStream = Files.newOutputStream(p);
            } else if (getURI().getScheme().equals(HDFS_SCHEME)) {
                org.apache.hadoop.fs.Path fsPath = new org.apache.hadoop.fs.Path(getURIString());
                org.apache.hadoop.fs.FileSystem fs = fsPath.getFileSystem(new org.apache.hadoop.conf.Configuration());
                return fs.create(fsPath);
            } else {
                outputStream = super.getOutputStream();
            }

            //TODO: should this wrap the stream in a .gz in a BlockCompressedOutputStream ?
            return outputStream;
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(getRawInputString(), e);
        }
    }

    @Override
    public void setTag(String tagName) {
        this.tagName = tagName;
    }

    @Override
    public String getTag() {
        return tagName;
    }

    @Override
    public void setTagAttributes(final Map<String, String> attributes) {
            this.tagAttributes = attributes;
    }

    @Override
    public Map<String, String> getTagAttributes() {
        return tagAttributes;
    }

    @Override
    public String toString() {
        if (getTagAttributes() != null) {
            return super.toString() + TaggedArgumentParser.getDisplayString("", this);
        } else {
            return super.toString();
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GATKPathSpecifier)) return false;
        if (!super.equals(o)) return false;

        GATKPathSpecifier GATKPathSpecifier = (GATKPathSpecifier) o;

        if (!Objects.equals(tagName, GATKPathSpecifier.tagName)) return false;
        return Objects.equals(getTagAttributes(), GATKPathSpecifier.getTagAttributes());
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = tagName == null ?
                result :
                31 * result + tagName.hashCode();
        result = getTagAttributes() == null ?
                result :
                31 * result + getTagAttributes().hashCode();
        return result;
    }
}
