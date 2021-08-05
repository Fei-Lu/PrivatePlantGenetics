/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.RNASeq;

import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.table.RowTable;
import gnu.trove.set.hash.TIntHashSet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
class Miscellaneous {
    
    public Miscellaneous () {
//        this.GFFToPGF();
//        this.barcodeTest();
        this.parseBarcodeTest();
    }

    public void parseBarcodeTest () {
        String barcodeFileS = "/Users/feilu/Documents/analysisL/pipelineTest/wega/parsingFastq/SampleInformation.txt";
        String fq1 = "/Users/feilu/Documents/analysisL/pipelineTest/wega/parsingFastq/20201216SIPASR2_R1.fq.gz";
        String fq2 = "/Users/feilu/Documents/analysisL/pipelineTest/wega/parsingFastq/20201216SIPASR2_R2.fq.gz";
        String subFqDirS = "/Users/feilu/Documents/analysisL/pipelineTest/wega/parsingFastq/subFq";
//        String barcodeFileS = "/data1/home/feilu/SampleInformation.txt";
//        String fq1 = "/data2/junxu/WEGAData/3/3_R1.fq.gz";
//        String fq2 = "/data2/junxu/WEGAData/3/3_R2.fq.gz";
//        String subFqDirS = "/data1/home/feilu/sub";

        new File (subFqDirS).mkdir();
        RowTable<String> t = new RowTable<>(barcodeFileS);
        HashMap<String, BufferedWriter> boMap1 = new HashMap();
        HashMap<String, BufferedWriter> boMap2 = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String outfileS1 = new File (subFqDirS, t.getCell(i,1)+"_r1.fq.gz").getAbsolutePath();
            String outfileS2 = new File (subFqDirS, t.getCell(i,1)+"_r2.fq.gz").getAbsolutePath();
            boMap1.put(t.getCell(i,6), IOUtils.getTextGzipWriter(outfileS1));
            boMap2.put(t.getCell(i,6), IOUtils.getTextGzipWriter(outfileS2));
        }
        Set<String> barcodeSet = boMap1.keySet();
        long startTime = System.nanoTime();
        String temp = null;String seq = null;
        try {
            BufferedReader br1 = IOUtils.getTextGzipReader(fq1);
            BufferedReader br2 = IOUtils.getTextGzipReader(fq2);

            String currentBarcode = null; String currentUMI = null;
            BufferedWriter tw1 = null;BufferedWriter tw2 = null;
            int cnt = 0;int cnt2 = 0;
            while((temp = br2.readLine())!=null){
                cnt2++;
                seq = br2.readLine();
                //*************************
                //barcode needs  to have at least 2 mismathches between each other, currently at 3 mismathces between any 8-bp barcodes
                //barcode can be redesigned to have 4-8 bp in length for even efficiency between barcdes
                //*************************
                currentBarcode = seq.substring(0, 8);
                tw2 = boMap2.get(currentBarcode);
                if (tw2 != null) {
                    tw2.write(temp);tw2.newLine();
                    tw2.write(seq);tw2.newLine();
                    tw2.write(br2.readLine());tw2.newLine();
                    tw2.write(br2.readLine());tw2.newLine();
                    tw1 = boMap1.get(currentBarcode);
                    tw1.write(br1.readLine());tw1.newLine();
                    tw1.write(br1.readLine());tw1.newLine();
                    tw1.write(br1.readLine());tw1.newLine();
                    tw1.write(br1.readLine());tw1.newLine();
                    cnt++;
                } else{
                    br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    br2.readLine();br2.readLine();
                }
                if (cnt2%1000000 == 0) System.out.println(cnt2);
            }
            StringBuilder sb = new StringBuilder();
            sb.append(cnt).append(" out of ").append(cnt2).append(", ").append(((float)(double)cnt/cnt2)).append(" of total reads were parsed from " ).append(" ").append(fq1);
            System.out.println(sb.toString());
            List<String> barcodeList = new ArrayList<>(barcodeSet);
            for (int i = 0; i < barcodeList.size(); i++) {
                tw1 = boMap1.get(barcodeList.get(i));
                tw2 = boMap2.get(barcodeList.get(i));
                tw1.flush();tw1.close();
                tw2.flush();tw2.close();
            }
            br1.close();
            br2.close();
            System.out.println(Benchmark.getTimeSpanMilliseconds(startTime));
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println(seq);
            System.exit(1);
        }
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
