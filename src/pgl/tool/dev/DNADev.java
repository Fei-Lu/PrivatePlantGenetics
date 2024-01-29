package pgl.tool.dev;

import pgl.infra.dna.FastaBit;
import pgl.infra.dna.FastaByte;
import pgl.infra.dna.allele.Allele;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;

import java.io.File;

public class DNADev {
    public DNADev () {
//        this.fasta();
    }

    public void fastq () {
        String infileS1 = "/Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/rawFastq/K16BJS0001_1.fq";
        String infileS2 = "/Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/rawFastq/K16BJS0001_2.fq";
    }

    public void fasta() {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa";
        String a = Benchmark.getMD5Checksum(infileS);
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/out.fa";
        FastaBit f = new FastaBit(infileS);
        f.writeFasta(outfileS, IOFileFormat.Text);
        String b = Benchmark.getMD5Checksum(outfileS);
        FastaByte fb = new FastaByte(infileS);
        fb.writeFasta(outfileS, IOFileFormat.Text);
        String c = Benchmark.getMD5Checksum(outfileS);
        System.out.println(a);
        System.out.println(b);
        System.out.println(c);
        new File(outfileS).delete();
    }
}
