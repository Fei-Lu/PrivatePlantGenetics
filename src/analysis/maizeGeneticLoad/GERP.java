/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maizeGeneticLoad;

import format.table.RowTable;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.CrossMapUtils;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author feilu
 */
class GERP {
    
    public GERP () {
        this.convertToV4GerpFile();
    }
    
    public void convertToV4GerpFile () {
        //this.mkAGPV3BED();
        this.crossMapConvert();
    }
    
    public void crossMapConvert () {
        String intDirS = "/Users/feilu/Documents/database/maize/crossMap/temp";
        String outDirS = "/Users/feilu/Documents/database/maize/crossMap/temp2";
        new File (outDirS).mkdir();
        File[] fs = new File(intDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "bed");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName()).getAbsolutePath();
            CrossMapUtils cm = new CrossMapUtils(f.getAbsolutePath(), outfileS);
            cm.setMaizeV3ToV4();
            cm.convert();
        });
    }
    
    public void mkAGPV3BED () {
        String infileS = "/Users/feilu/Documents/database/maize/infoFile/ChrLenCentPosi_agpV3.txt";
        String outDirS = "/Users/feilu/Documents/database/maize/crossMap/temp";
        RowTable<String> t = new RowTable(infileS);
        List<Integer> chrList = new ArrayList();
        int[] chrLength = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            String s = t.getCell(i, 0);
            chrList.add(Integer.parseInt(s));
            chrLength[i] = Integer.parseInt(t.getCell(i, 1));
        }
        chrList.parallelStream().forEach(chr -> {
            String outfileS = new File (outDirS, "chr"+PStringUtils.getNDigitNumber(3, chr)+".bed").getAbsolutePath();
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < chrLength[chr-1]; i++) {
                    sb = new StringBuilder();
                    sb.append(chr).append("\t").append(i+1).append("\t").append(i+1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
}
