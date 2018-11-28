package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

public class ExampleReferenceWalkerIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testExampleReferenceWalker(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg19MiniReference))
                .addArgument("L", "1:1000-2000");
        runCommandLine(args);
    }


    @Test
    public void testExampleReferenceWalkerFull(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference));
        runCommandLine(args);
    }
}