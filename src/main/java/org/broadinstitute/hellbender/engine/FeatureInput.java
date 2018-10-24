package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.Serializable;
import java.net.URI;
import java.util.*;

/**
 * Class to represent a Feature-containing input file. Tools should declare @Argument-annotated fields of
 * this type (or Collections of this type), and the Feature management system will automatically discover
 * them at runtime (provided that they are declared in the tool itself, a superclass of the tool, or an
 * ArgumentCollection of the tool).
 *
 * DO NOT ATTEMPT TO INSTANTIATE THIS CLASS DIRECTLY! FeatureInputs must be instantiated by the argument-parsing
 * system only in order to be recognized by the Feature management system. This is why the constructor is
 * marked as protected.
 *
 * If you still want to instantiate this class directly, you will have to call {@link GATKTool#addFeatureInputsAfterInitialization(String, String, Class, int)}
 *  in order to register the FeatureInput with the engine.
 *
 * FeatureInputs can be assigned logical names on the command line using the syntax:
 *
 *     --argument_name:logical_name feature_file
 *
 * These logical names can then be retrieved by the tool at runtime via {@link #getName}
 *
 * Furthermore, a list of comma-separated key=value pairs may be provided as follows:
 *
 *     --argument_name:logical_name,key1=value1,key2=value2 feature_file
 *
 * the string value provided for a given key can be retrieved via {@link #getAttribute(String)}. Keys must be unique.
 *
 * @param <T> the type of Feature that this FeatureInput file contains (eg., VariantContext, BEDFeature, etc.)
 */
public final class FeatureInput<T extends Feature> extends GATKInputPath implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * File containing Features as specified by the user on the command line
     */
    private final String rawInputSpecifier;

    /**
     * Cache the codec for this feature input the first time we discover it, so we only do it once
     */
    private transient Class<FeatureCodec<T, ?>> featureCodecClass;

    /**
     * Delimiter between the logical name and the file name in the --argument_name logical_name:feature_file syntax
     */
    public static final String FEATURE_ARGUMENT_TAG_DELIMITER = ":";

    /**
     * Construct a FeatureInput from a String argument value.
     *
     * Only meant to be called by the argument parsing system, and therefore marked as package-visible --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     *
     * Note: cannot delegate to another constructor because Java only allows a call to "this" on the first line of a constructor.
     *
     * @param rawArgumentValue String of the form "logical_name:feature_file" or "feature_file"
     */
    FeatureInput(final String rawArgumentValue) {
        super(rawArgumentValue);
        Utils.nonNull(rawArgumentValue, "rawArgumentValue");
        this.rawInputSpecifier = rawArgumentValue;
        // TODO restore this
//        if (getTag() == null) {
//            setTag(rawArgumentValue);
//        }
    }

    /**
     * Construct a FeatureInput from a path and a name
     *
     * This constructor is meant to be called only by the engine and test classes,
     * which is why it has package access.
     */
    FeatureInput(final String featurePath, final String name) {
        this(featurePath, name, Collections.emptyMap());
    }

    /**
     * Construct a FeatureInput from raw components: name, key value pairs and the file.
     *
     * This constructor is meant to be called by the engine and test classes --
     * FeatureInputs constructed some other way will not be recognized by the engine.
     */
    @VisibleForTesting
    public FeatureInput(final String rawInputSpecifier, final String name, final Map<String, String> keyValueMap) {
        super(rawInputSpecifier);

        Utils.nonNull(name, "name");
        Utils.nonNull(keyValueMap, "kevValueMap");
        Utils.nonNull(rawInputSpecifier, "feature-file");

        setTag(name);
        setTagAttributes(keyValueMap);
        this.rawInputSpecifier = rawInputSpecifier;
    }

    /**
     * Remember the FeatureCodec class for this input the first time it is discovered so we can bypass dynamic codec
     * discovery when multiple FeatureDataSources are created for the same input.
     */
    public void setFeatureCodecClass(final Class<FeatureCodec<T, ?>> featureCodecClass) {
        this.featureCodecClass = featureCodecClass;
    }

    /**
     * @return The previously established FeatureCodec class to use for this input, if any. May return {@code null}.
     */
    public Class<FeatureCodec<T, ?>> getFeatureCodecClass() {
        return this.featureCodecClass;
    }

    /**
     * creates a name from the given filePath by finding the absolute path of the given input
     */
    //TODO: this is redundant with BucketUtils.makeFilePathAbsolute(String)
    private static String makeIntoAbsolutePath(final String filePath){
        // TODO: update these to use GatkUri
        if(IOUtils.isGenomicsDBPath(filePath)){
            return IOUtils.getAbsolutePathWithGenomicsDBURIScheme(filePath);
        } else if (URI.create(filePath).getScheme() != null) {
            return IOUtils.getPath(filePath).toAbsolutePath().toUri().toString();
        } else {
            return new File(filePath).getAbsolutePath();
        }
    }

    /**
     * Gets the value for the given key associated with this Feature source or {@code null}
     * if no value is associated with a given key.
     * @throws IllegalArgumentException if the key is {@code null}.
     */
    public String getAttribute(final String key) {
        Utils.nonNull(key);
        return super.getTagAttributes().get(key);
    }

    /**
     * Gets the logical name of this Feature source. This will be a user-provided value if the
     * --argument_name logical_name:feature_file was used on the command line, otherwise it will
     * default to the absolute path of the backing file
     *
     * @return logical name of this source of Features
     */
    public String getName() {
        return super.getTag() != null ? super.getTag() : getURIString();
    }

    /**
     *
     * @return true if the value for name does not match the default, indicating it was a user supplied name (i.e. foo:file.vcf)
     */
    public boolean hasUserSuppliedName() { return getTag() != null;
    }

    /**
     * Gets the file backing this source of Features
     *
     * @return file backing this source of Features
     */
    //TODO: consumers of this method assume this is a File
    public String getFeaturePath() {
        return super.getURI().getPath();
    }

    /**
     * FeatureInputs will be hashed by the engine, so make an effort to produce a reasonable hash code
     *
     * @return hash code for this FeatureInput (combination of hash code of the name and file)
     */
    @Override
    public int hashCode() {
        return super.hashCode() + 31 * rawInputSpecifier.hashCode();
    }

    /**
     * Test FeatureInputs for equality
     *
     * @param other object to compare this FeatureInput with
     * @return true if this FeatureInput equals other, otherwise false
     */
    @Override
    public boolean equals(final Object other) {
        if (! (other instanceof FeatureInput)) {
            return false;
        }

        final FeatureInput<?> otherFeature = (FeatureInput<?>)other;
        if (!Objects.equals(getTag(), otherFeature.getTag())) {
            return false;
        }

        return rawInputSpecifier.equals(otherFeature.rawInputSpecifier);
    }

    /**
     * Returns a String representation of this FeatureInput. Will be the absolute path to
     * the featureFile if we have no logical name, or a String of the form
     * "logical_name:absolute_path_to_featureFile" if we do have a logical name.
     *
     * @return String representation of this FeatureInput
     */
    @Override
    public String toString() {
        //TODO: do we need this call to makeIntoAbsolutePath ?
//        final String featureFilePath = makeIntoAbsolutePath(rawInputSpecifier);
//        return getTag() != null && featureFilePath.equals(getTag()) ?
//                rawInputSpecifier :
//                String.format("%s%s%s", getTag(), FEATURE_ARGUMENT_TAG_DELIMITER, featureFilePath);
        return rawInputSpecifier;
    }
}
