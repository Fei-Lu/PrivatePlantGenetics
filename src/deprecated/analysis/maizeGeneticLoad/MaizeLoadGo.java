/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deprecated.analysis.maizeGeneticLoad;

import pgl.format.table.RowTable;
import java.util.List;

/**
 *
 * @author feilu
 */
public class MaizeLoadGo {
    
    public MaizeLoadGo () {
        this.analysisPipe();
    }
    
    public void analysisPipe () {
        //this.annotateVariants();
        //this.processGERP();
        this.processSIFT();
        //test();
    }
    
    public void test () {
        String inputFileS = "/Users/feilu/Desktop/hmp321Info_chr010.txt";
        String inputBedFileS = "/Users/feilu/Documents/database/maize/crossMap/test/test.bed";
        RowTable<String> t = new RowTable(inputFileS);
        int[] chr = new int[t.getRowNumber()];
        int[] pos = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            chr[i] = Integer.parseInt(t.getCell(i, 0));
            pos[i] = Integer.parseInt(t.getCell(i, 1));
        }
        CrossMapUtils cm = new CrossMapUtils(chr, pos, inputBedFileS);
        cm.setMaizeV3ToV4();
        cm.convert();
        List<int[]> l  = cm.getConvertedCoordinate();
        int[] newPos = l.get(1);
        int cnt = 0;
        for (int i = 0; i < newPos.length; i++) {
            if (newPos[i] == -1) continue;
            cnt++;
        }
        System.out.println(cnt);
    }
    
    public void annotateVariants () {
        new VariantsAnnotation ();
    }
    
    public void processGERP () {
        new GERP();
    }
    
    public void processSIFT () {
        new SIFT ();
    }
    
    public static void main (String[] args) {
        new MaizeLoadGo ();
    }
    
}
