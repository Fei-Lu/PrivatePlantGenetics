/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.RNASeq;

import format.dna.FastaBit;
import format.genomeAnnotation.GeneFeature;

/**
 *
 * @author feilu
 */
class Miscellaneous {
    
    public Miscellaneous () {
        this.GFFToPGF();
    }
    
    
    public void GFFToPGF () {
        String inputGFF = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3";
        String outputPGF = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature();
        gf.readFromWheatGFF(inputGFF);
        gf.writeFile(outputPGF);
    }
    
}