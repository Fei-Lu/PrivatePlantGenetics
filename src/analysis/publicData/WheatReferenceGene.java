/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.publicData;

import format.genomeAnnotation.GeneFeature;
import format.table.RowTable;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import utils.Benchmark;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class WheatReferenceGene {
    
    public WheatReferenceGene () {
        //this.GFFToPGF();
        //this.mkExpressionByTissue();
        this.mkExpressedGene();
    }
    
    public void mkExpressedGene () {
        String pgfFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String expressionInFileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/download/Development_tpm.tsv.gz";
        String header = "Gene\tOverlapped_gene\tLongest_overlapped_gene\tLongest_transcript\tTPM_mean\tTPM_sd\tTPM_rsd";
        
    }
    
    public void mkExpressionByTissue () {
        String sampleInFileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/download/Table_S1.txt";
        String expressionInFileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/download/Development_tpm.tsv.gz";
        String outfileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/expressionByTissue.txt";
        RowTable<String> t = new RowTable(sampleInFileS);
        HashSet<String> s = new HashSet(t.getColumn(9));
        String[] tissues = s.toArray(new String[s.size()]);
        Arrays.sort(tissues);
        List<String>[] sampleLists = new List[tissues.length];
        for (int i = 0; i < sampleLists.length; i++) sampleLists[i] = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(tissues, t.getCell(i, 9));
            sampleLists[index].add(t.getCell(i, 0));
        }
        int[][] sampleIndex = new int[sampleLists.length][];
        for (int i = 0; i < sampleLists.length; i++) {
            sampleIndex[i] = new int[sampleLists[i].size()];
            Collections.sort(sampleLists[i]);
        }
        t = new RowTable(expressionInFileS);
        List<String> header = t.getHeader();
        for (int i = 0; i < header.size(); i++) {
            for (int j = 0; j < tissues.length; j++)  {
                int index = Collections.binarySearch(sampleLists[j], header.get(i));
                if (index < 0) continue;
                sampleIndex[j][index] = i;
            }
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("Transcripts");
            for (int i = 0; i < tissues.length; i++) {
                sb.append("\t").append(tissues[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append(t.getCell(i, 0));
                for (int j = 0; j < tissues.length; j++) {
                    double e = 0;
                    for (int k = 0; k < sampleLists[j].size(); k++) {
                        e+=Double.parseDouble(t.getCell(i, sampleIndex[j][k]));
                    }
                    e = e/sampleIndex[j].length;
                    sb.append("\t").append((float)e);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void GFFToPGF () {
        String inputGFF = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3";
        String outputPGF = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature();
        gf.readFromWheatGFF(inputGFF);
        gf.writeFile(outputPGF);
    }
    
}
