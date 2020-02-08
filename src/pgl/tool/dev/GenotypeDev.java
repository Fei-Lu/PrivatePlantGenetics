package pgl.tool.dev;

import pgl.infra.dna.allele.AlleleEncoder;

public class GenotypeDev {

    public GenotypeDev () {
        this.testAlleles();
    }

    public void testAlleles () {
        String genotype = "A.";
        byte b1 = AlleleEncoder.getAlleleByteFromChar(genotype.charAt(0));
        byte b2 = AlleleEncoder.getAlleleByteFromChar(genotype.charAt(1));
        byte geno = AlleleEncoder.getGenotypeByte(b1, b2);
        System.out.println(geno);
        System.out.println(AlleleEncoder.getAlleleChar1FromGenotypeByte(geno));
        System.out.println(AlleleEncoder.getAlleleChar2FromGenotypeByte(geno));
    }
}
