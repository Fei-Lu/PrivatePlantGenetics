package pgl.tool.dev;

import pgl.infra.dna.allele.Allele;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.genotype.GenoIOFormat;
import pgl.infra.dna.genotype.GenotypeBit;
import pgl.infra.dna.genotype.GenotypeExport;
import pgl.infra.dna.genotype.GenotypeTable;
import pgl.infra.utils.Benchmark;

import java.nio.ByteBuffer;
import java.util.BitSet;

public class GenotypeDev {

    public GenotypeDev () {
        //this.testAlleles();
        //this.ioTest();
        this.dxyTest();
    }

    public void dxyTest () {
        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/out.bin.gz";
        GenotypeTable gt  = new GenotypeBit(vcfFileS, GenoIOFormat.Binary_GZ);
        System.out.println(gt.getTaxaNumber());
        System.out.println(gt.getSiteNumber());
        long start = System.nanoTime();
        double[][] matrix = gt.getDxyMatrixFast10K();
        System.out.println(Benchmark.getTimeSpanMilliseconds(start));
    }

    public void ioTest () {
        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon_vmap2.1.vcf.gz";
        GenotypeTable gt  = new GenotypeBit(vcfFileS, GenoIOFormat.VCF_GZ);

        String vcfOutGZFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/out.vcf.gz";
        GenotypeExport.output(gt, vcfOutGZFileS, GenoIOFormat.VCF_GZ);

        String binOutFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/out.bin.gz";
        GenotypeExport.output(gt, binOutFileS, GenoIOFormat.Binary_GZ);

//        String vcfOutFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/out.vcf";
//        gt  = new GenotypeBit(binOutFileS, GenoIOFormat.Binary_GZ);
//        GenotypeExport.output(gt, vcfOutFileS, GenoIOFormat.VCF);

    }

    public void testAlleles () {
        String genotype = "A.";
        byte b1 = AlleleEncoder.getAlleleByteFromBase(genotype.charAt(0));
        byte b2 = AlleleEncoder.getAlleleByteFromBase(genotype.charAt(1));
        byte geno = AlleleEncoder.getGenotypeByte(b1, b2);
        System.out.println(geno);
        System.out.println(AlleleEncoder.getAlleleBase1FromGenotypeByte(geno));
        System.out.println(AlleleEncoder.getAlleleBase2FromGenotypeByte(geno));

        Allele al = new Allele ('A');
        al.setAlleleType(AlleleType.Minor);
        al.setAlleleType(AlleleType.Reference);
        if (al.isAlleleTypeOf(AlleleType.Minor)) {
            System.out.println("yes");
        }
        System.out.println(Integer.toBinaryString(al.getAlleleFeature()));
    }
}
