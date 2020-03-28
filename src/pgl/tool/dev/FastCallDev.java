package pgl.tool.dev;

import pgl.infra.dna.*;
import pgl.infra.utils.IOFileFormat;

public class FastCallDev {

    public FastCallDev () {
        this.mkSubReference();
    }

    public void mkSubReference () {
        int subGenomeSize = 1_000_000;
        String originalRef = "/Users/feilu/Documents/database/wheat/reference/v1.0/byChr/chr001.fa.gz";
        String targetRef = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa.gz";
        FastaBit fa = new FastaBit(originalRef);
        FastaRecordBit fb = fa.getFastaRecordBit(0, 0, subGenomeSize);
        fa = new FastaBit(fb);
        fa.writeFasta(targetRef, IOFileFormat.TextGzip);
    }



}

