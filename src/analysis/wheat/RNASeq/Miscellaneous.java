/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.RNASeq;

import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.table.RowTable;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedWriter;
import java.util.Arrays;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
class Miscellaneous {
    
    public Miscellaneous () {
//        this.GFFToPGF();
        this.barcodeTest();
    }
    
    public void barcodeTest () {
        String barcodeFileS = "/Users/feilu/Downloads/barcode.txt";
        String outFileS = "/Users/feilu/Downloads/out.txt";
        RowTable<String> t = new RowTable (barcodeFileS);
        TIntHashSet tL = new TIntHashSet();
        TIntHashSet tS = new TIntHashSet();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outFileS);
            bw.write("IDstart\tIDEnd\tDiff");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                for (int j = 0; j < t.getRowNumber(); j++) {
                    if (i == j) continue;
                    StringBuilder sb = new StringBuilder();
                    byte[] s1 = t.getCell(i, 2).getBytes();
                    byte[] s2 = t.getCell(j, 2).getBytes();
                    int cnt = 0;
                    for (int k = 0; k < s1.length; k++) {
                        if (s1[k] != s2[k]) cnt++;
                    }
                    if (cnt == 0) {
                        
                    }
                    else if (cnt > 2) {
                        tL.add(i);
                        tL.add(j);
                    }
                    else {
                        tS.add(i);
                        tS.add(j);
                    }
                    sb.append(i).append("\t").append(j).append("\t").append(cnt);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int a = 3;
        int[] tsA = tS.toArray();
        for (int i = 0; i < tS.size(); i++) {
            if (tL.contains(tsA[i])) {
                tL.remove(tsA[i]);
            }
        }
        int[] tsL = tL.toArray();
        Arrays.sort(tsL);
                
        for (int i = 0; i < tsL.length; i++) {
            System.out.println(tsL[i]);
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
