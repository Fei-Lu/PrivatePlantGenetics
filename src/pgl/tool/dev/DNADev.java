package pgl.tool.dev;

import pgl.infra.dna.FastaBit;
import pgl.infra.dna.allele.Allele;

public class DNADev {
    public DNADev () {
        this.testFasta();
    }

    public void testFasta () {
        String infileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/byChr/chr044.fa.gz";
        FastaBit fb = new FastaBit(infileS);

        Allele al = new Allele ('T');
    }
}
