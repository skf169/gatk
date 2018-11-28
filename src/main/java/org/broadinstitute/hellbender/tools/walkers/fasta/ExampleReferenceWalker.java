package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceWalker;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.util.HashMap;
import java.util.Map;

@CommandLineProgramProperties(
  summary = "Count the number of bases in a reference file",
  oneLineSummary = "Count the number of bases that appear in a reference file",
  programGroup = ReferenceProgramGroup.class
)
public class ExampleReferenceWalker extends ReferenceWalker {
    private Map<String, Long> contigBaseCount = new HashMap<>();

    @Override
    public void apply(ReferenceContext referenceContext, ReadsContext read, FeatureContext featureContext) {
      contigBaseCount.merge(referenceContext.getInterval().getContig(), 1L, (old, newValue) -> old + newValue);
    }

    @Override
    public Object onTraversalSuccess(){
      contigBaseCount.forEach((key, value) -> System.out.println(key + " : " + value));
      final long totalLength = contigBaseCount.values().stream().mapToLong(Long::longValue).sum();
      System.out.println(totalLength);
      return totalLength;
    }
}
