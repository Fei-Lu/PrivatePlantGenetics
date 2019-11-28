/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.tablesaw.TablesawUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;
import utils.IOUtils;
import utils.PStringUtils;
/**
 *
 * @author feilu
 */
public class DeleteriousDB {
    
    public DeleteriousDB () {
       //this.extractInfoFromVMap2();
       this.mkGenicAnnotation();
       this.addSift();
    }
    
    public void addSift () {
        String siftDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/sift/output/";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation";
        File[] fs = new File(siftDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".xls");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            //f = new File("/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/sift/output/chr001.subgenome.maf0.01byPop.SNP_SIFTannotations.xls");
            String dbFileS = f.getName().split("\\.")[0]+"_SNP_anno.txt";
            dbFileS = new File (dirS, dbFileS).getAbsolutePath();
            Table td = TablesawUtils.readTsv(dbFileS);
            Table ts = TablesawUtils.readTsv(f.getAbsolutePath());
            ts.sortAscendingOn("POS");
            String[] region = new String[td.rowCount()];
            String[] type = new String[td.rowCount()];
            String[] sift = new String[td.rowCount()];
            int pos = 0;
            String trans = null;
            String alt = null;
            int cnt = 0;
            for (int i = 0; i < td.rowCount(); i++) {
                pos = td.intColumn("Pos").getInt(i);
                trans = td.getString(i, "Transcript");
                alt = td.getString(i, "Alt");
                Table result = ts.where(ts.intColumn("POS").isEqualTo(pos).and(ts.stringColumn("TRANSCRIPT_ID").isEqualTo(trans)));
                if (result.isEmpty()) {
                    region[i] = "NA";
                    type[i] = "NONCODING";
                    sift[i] = "NA";
                }
                else {
                    Table result2 = result.where(result.stringColumn("ALT_ALLELE").isEqualTo(alt));
                    if (result2.isEmpty()) {
                        
                    }
                    else {
                        region[i] = result2.getString(0, "REGION");
                        type[i] = result2.getString(0, "VARIANT_TYPE");
                        sift[i] = result2.getString(0, "SIFT_SCORE");
                        if (sift[i].isEmpty()) sift[i] = "NA";
                    }
                }
                cnt++;
                if (cnt%10000 == 0) System.out.println(cnt+"\t"+dbFileS);
            }
            StringColumn regionC = StringColumn.create("Region", region);
            StringColumn typeC = StringColumn.create("Variant_type", type);
            StringColumn siftC = StringColumn.create("SIFT_score", sift);
            int columnIndex = td.columnCount();
            td.insertColumn(columnIndex, regionC);
            td.insertColumn(++columnIndex, typeC);
            td.insertColumn(++columnIndex, siftC);
            td = td.where(td.stringColumn("Region").isNotEqualTo("NA"));
            TablesawUtils.writeTsv(td, dbFileS);
            System.out.println(f.getAbsolutePath()+ " completed");
        });
        
    }
    
    public void mkGenicAnnotation () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPByChr/";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = f.getName().split("_")[0]+"_SNP_anno.txt";
            outfileS = new File (outDirS, outfileS).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
    
    public void extractInfoFromVMap2 () {
        int subLength = 150;
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPByChr/";
        String vmapDirS = "/Volumes/Fei_HDD_Mac/VMap2.1/";
        File[] fs  = new File(vmapDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        String geneHCFileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        Table t = TablesawUtils.readTsv(geneHCFileS);
        System.out.println(t.structure());
        t.sortAscendingOn("Chr", "TranStart");
        IntColumn chrColumn = t.intColumn("chr");
        int chrNum = chrColumn.countUnique();
        TIntList[] startLists = new TIntList[chrNum];
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }
        for (int i = 0; i < t.rowCount(); i++) {
            startLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 3)));
            endLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 4)));
            tranLists[Integer.parseInt(t.getString(i, 2))-1].add(t.getString(i, 1));
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst(".vcf.gz", "genicSNP.txt.gz")).getAbsolutePath();
            int[] dc = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(dc);
            StringBuilder sb = new StringBuilder();
            if (Arrays.binarySearch(dc, chrIndex+1) < 0) {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB\tTranscript");
            }
            else {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D\tTranscript");
            }
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()).startsWith("#")) {}
                
                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) {
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    posIndex = startLists[chrIndex].binarySearch(currentPos);
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    sb.append("\t").append(l.get(4)).append("\t");
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""),",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    }
                    else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]).append("\t").append(ll.get(7).split("=")[1]).append("\t").append(ll.get(8).split("=")[1]);
                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
}
