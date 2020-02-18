package pgl.tool.dev;

import pgl.infra.dna.allele.Allele;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.genotype.GenoIOFormat;
import pgl.infra.dna.genotype.GenotypeBit;
import pgl.infra.dna.genotype.GenotypeTable;

public class GenotypeDev {

    public GenotypeDev () {
        //this.testAlleles();
        this.readVCF();
    }

    public void readVCF () {
        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/chr001_exon_vmap2.1.vcf.gz";
        GenotypeTable gt  = new GenotypeBit(vcfFileS, GenoIOFormat.VCF);
        System.out.println(gt.getMinorAlleleFrequency(1));
        gt.sortByTaxa();
        System.out.println(gt.getMinorAlleleFrequency(1));
        

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
