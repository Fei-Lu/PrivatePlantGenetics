package pgl.tool.dev;

import pgl.infra.genomeAnnotation.GeneFeature;

public class GenomeAnnotationDev {
    public GenomeAnnotationDev () {
        this.geneFeature();
    }

    public void geneFeature () {
        String infileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature(infileS);
        int index = gf.getGeneIndex(0,	1147);
        System.out.println(index);
    }
}
