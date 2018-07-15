/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maizeRNASeq;

import format.genomeAnnotation.GeneFeature;
import format.range.Range;
import format.range.Ranges;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author feilu
 */
public class RNASeqMiscellaneous {
    
    public RNASeqMiscellaneous () {
        this.testOverlapGeneModel();
    }
    
    public void testOverlapGeneModel () {
        String infileS = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.pgf";
        GeneFeature gf = new GeneFeature(infileS);
        List<Range> rList = new ArrayList<>();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (gf.getGeneChromosome(i) != 1) continue;
            if (gf.getGeneStrand(i) != 1) continue;
            Range r = new Range (gf.getGeneChromosome(i), gf.getGeneStart(i), gf.getGeneEnd(i));
            rList.add(r);
        }
        Ranges rs = new Ranges(rList);
        System.out.println(rs.getRangeNumber());
        Ranges nrs = rs.getNonOverlapRanges();
        System.out.println(nrs.getRangeNumber());
    }
}
