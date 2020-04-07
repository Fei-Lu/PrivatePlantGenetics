package pgl.tool.dev;

import pgl.infra.dna.allele.Allele;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.genotype.*;
import pgl.infra.dna.genotype.summary.SumTaxaDivergence;
import pgl.infra.utils.IOFileFormat;

public class GenotypeDev {

    public GenotypeDev () {
//        this.testAlleles();
//        this.ioTestGenotypeBit();
        this.ioTestGenotypeGrid();
//        this.dxyTest1();
//        this.dxyTest2();
    }

    public void dxyTest2 () {
        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/source/chr001_exon_vmap2.1.vcf.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/taxaDivergence.txt";
        GenotypeTable gt  = new GenotypeGrid(vcfFileS, GenoIOFormat.VCF_GZ);
        SumTaxaDivergence std = new SumTaxaDivergence(gt);
        std.writeSummary(outfileS, IOFileFormat.Text);
        gt.getDxyMatrix();
    }

    public void dxyTest1 () {
//        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon.bin.gz";
//        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/taxaDivergence.txt";
//        GenotypeTable gt  = new GenotypeBit(vcfFileS, GenoIOFormat.Binary_GZ);
//
//        System.out.println(gt.getTaxaNumber());
//        System.out.println(gt.getSiteNumber());
//
//        SumTaxaDivergence std = new SumTaxaDivergence(gt);
//        std.writeSummary(outfileS, IOFileFormat.Text);

    }

    public void ioTestGenotypeGrid() {

        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/source/chr001_exon_vmap2.1.vcf.gz";
        GenotypeTable gt  = new GenotypeGrid(vcfFileS, GenoIOFormat.VCF_GZ);

        String vcfOutGZFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon.vcf";
        GenotypeExport.output(gt, vcfOutGZFileS, GenoIOFormat.VCF);

        String binOutFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon.bin.gz";
        GenotypeExport.output(gt, binOutFileS, GenoIOFormat.Binary_GZ);

        gt  = new GenotypeGrid(binOutFileS, GenoIOFormat.Binary_GZ);
        String vcffromBinFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon_fromBin.vcf";
        GenotypeExport.output(gt, vcffromBinFileS, GenoIOFormat.VCF);

        gt  = new GenotypeGrid(vcfOutGZFileS, GenoIOFormat.VCF);
        String vcffromVcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon_fromVCF.vcf";
        GenotypeExport.output(gt, vcffromVcfFileS, GenoIOFormat.VCF);
    }

    public void ioTestGenotypeBit() {
//        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/source/chr001_exon_vmap2.1.vcf.gz";
//        GenotypeTable gt  = new GenotypeBit(vcfFileS, GenoIOFormat.VCF_GZ);
//
//        String vcfOutGZFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon.vcf";
//        GenotypeExport.output(gt, vcfOutGZFileS, GenoIOFormat.VCF);
//
//        String binOutFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon.bin.gz";
//        GenotypeExport.output(gt, binOutFileS, GenoIOFormat.Binary_GZ);
//
//        gt  = new GenotypeBit(binOutFileS, GenoIOFormat.Binary_GZ);
//        String vcffromBinFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon_fromBin.vcf";
//        GenotypeExport.output(gt, vcffromBinFileS, GenoIOFormat.VCF);
//
//        gt  = new GenotypeBit(vcfOutGZFileS, GenoIOFormat.VCF);
//        String vcffromVcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon_fromVCF.vcf";
//        GenotypeExport.output(gt, vcffromVcfFileS, GenoIOFormat.VCF);
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
