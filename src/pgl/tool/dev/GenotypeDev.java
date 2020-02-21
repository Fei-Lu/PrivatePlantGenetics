package pgl.tool.dev;

import pgl.infra.dna.allele.Allele;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.genotype.GenoIOFormat;
import pgl.infra.dna.genotype.GenotypeBit;
import pgl.infra.dna.genotype.GenotypeExport;
import pgl.infra.dna.genotype.GenotypeTable;

import java.nio.ByteBuffer;
import java.util.BitSet;

public class GenotypeDev {

    public GenotypeDev () {
        //this.testAlleles();
        this.readVCF();
        BitSet bs = new BitSet(65);
        bs.set(65);

        long[] bl = bs.toLongArray();
        byte[] ba = bs.toByteArray();
        int v = bs.size();
        int a = 3;
        
    }

    public void readVCF () {
        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon_vmap2.1.vcf.gz";
//        String vcfFileS = "/Volumes/Fei_HDD_Mac/VMap2.1/chr001_vmap2.1.vcf.gz";
        GenotypeTable gt  = new GenotypeBit(vcfFileS, GenoIOFormat.VCF_GZ);

        String vcfOutFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/out.vcf";
        GenotypeExport.output(gt, vcfOutFileS, GenoIOFormat.VCF);

        String binOutFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/out.bin.gz";
        GenotypeExport.output(gt, binOutFileS, GenoIOFormat.Binary_GZ);

        gt  = new GenotypeBit(binOutFileS, GenoIOFormat.Binary_GZ);
        GenotypeExport.output(gt, vcfOutFileS, GenoIOFormat.VCF);

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
