package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;

public abstract class ReferenceWalker extends GATKTool {

  @Override
  public String getProgressMeterRecordLabel() { return "bases"; }

  @Override
  public boolean requiresReference() { return true; }

  /**
   * Initialize data sources for traversal.
   *
   * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
   */
  @Override
  protected final void onStartup() {
    super.onStartup();
  }


  /**
   * Implementation of read-based traversal.
   * Subclasses can override to provide own behavior but default implementation should be suitable for most uses.
   *
   * The default implementation creates filters using {@link #makeReadFilter} and transformers using
   * {@link #makePreReadFilterTransformer()} {@link #makePostReadFilterTransformer()} and then iterates over all reads, applies
   * the pre-filter transformer, the filter, then the post-filter transformer and hands the resulting reads to the {@link #apply}
   * function of the walker (along with additional contextual information, if present, such as reference bases).
   */
  @Override
  public void traverse() {
    // Process each read in the input stream.
    // Supply reference bases spanning each read, if a reference is available.
    final CountingReadFilter readFilter = makeReadFilter();



    for(SimpleInterval locus : getIntervalIterator()){
        final ReferenceContext referenceContext = new ReferenceContext(reference, locus, getWindowLeadingBases(), getWindowTrailingBases());
        apply(referenceContext,
              new ReadsContext(reads, referenceContext.getWindow(), readFilter), // Will create an empty ReadsContext if reads or readInterval == null
              new FeatureContext(features, referenceContext.getWindow()));   // Will create an empty FeatureContext if features or readInterval == null


        progressMeter.update(referenceContext.getInterval());
      };

  }

  protected int getWindowTrailingBases() {
    return 0;
  }

  protected int getWindowLeadingBases(){
    return 0;
  }

  private Iterable<SimpleInterval> getIntervalIterator(){
    return () -> new IntervalLocusIterator(getTraversalIntervals().iterator());
  }

  /**
   * Process an individual read (with optional contextual information). Must be implemented by tool authors.
   * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
   * as possible.
   *
   * TODO: Determine whether and to what degree the GATK engine should provide a reduce operation
   * TODO: to complement this operation. At a minimum, we should make apply() return a value to
   * TODO: discourage statefulness in walkers, but how this value should be handled is TBD.
   * @param read current read
   * @param referenceContext Reference bases spanning the current read. Will be an empty, but non-null, context object
   *                         if there is no backing source of reference data (in which case all queries on it will return
   *                         an empty array/iterator). Can request extra bases of context around the current read's interval
   *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
   * @param featureContext Features spanning the current read. Will be an empty, but non-null, context object
   *                       if there is no backing source of Feature data (in which case all queries on it will return an
   *                       empty List).
   */
  public abstract void apply(ReferenceContext referenceContext, ReadsContext read, FeatureContext featureContext );

  /**
   * Shutdown data sources.
   *
   * Marked final so that tool authors don't override it. Tool authors should override onTraversalSuccess() instead.
   */
  @Override
  protected final void onShutdown() {
    // Overridden only to make final so that concrete tool implementations don't override
    super.onShutdown();
  }
}
