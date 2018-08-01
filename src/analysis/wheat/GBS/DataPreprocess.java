/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.GBS;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class DataPreprocess {
    
    public DataPreprocess () {
        //this.parseFastqByIndex();
        this.mkNuoheLibFastqMapFileS();
    }
    
    public void mkNuoheLibFastqMapFileS () {
        String fqDirS = "/Users/feilu/Documents/analysisL/production/wheatGBS/source/nuohe";
        String barcodeFileS = "/Users/feilu/Documents/analysisL/production/wheatGBS/source/barcode_wheat_201807_modified.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/production/wheatGBS/source/libraryFastqMap_nuohe.txt";
        File[] subFiles = IOUtils.listRecursiveFiles(new File(fqDirS));
        File[] fqs = IOUtils.listFilesEndsWith(subFiles, ".gz");
        RowTable<String> t = new RowTable(barcodeFileS);
        List<String> l = t.getColumn(0);
        Set<String> lSet = new HashSet<>(l);
        String[] libs = lSet.toArray(new String[lSet.size()]);
        Arrays.sort(libs);
        String[] r1s = new String[libs.length];
        String[] r2s = new String[libs.length];
        for (int i = 0; i < fqs.length; i++) {
            String name = fqs[i].getName();
            String[] temp = name.split("_");
            String libName = temp[0].replaceFirst("s", "");
            int index = Arrays.binarySearch(libs, libName);
            if (index < 0) {
                continue;
            }
            if (temp[2].startsWith("1")) r1s[index] = fqs[i].getAbsolutePath();
            else if (temp[2].startsWith("2")) r2s[index] = fqs[i].getAbsolutePath();
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Library-name	R1Path	R2Path");
            bw.newLine();
            for (int i = 0; i < libs.length; i++) {
                StringBuilder sb = new StringBuilder();
                if (r1s[i] == null) r1s[i] = "NA";
                if (r2s[i] == null) r2s[i] = "NA";
                sb.append(libs[i]).append("\t").append(r1s[i]).append("\t").append(r2s[i]);
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
    
    public void parseFastqByIndex () {
        String infileS = "/Users/feilu/Documents/analysisL/production/wheatGBS/source/annuo/cleandata/20180705mix1/20180705mix1_R1.fq.gz";
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            for (int i = 0; i < 200; i++) {
                //if (i < 200) continue;
                System.out.println(br.readLine());
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
