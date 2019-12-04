/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import format.table.RowTable;
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
       //this.mkGenicAnnotation();
       //this.addSift2();
       //this.addAncestral();
       //this.addDAF();
       this.addGerp();
    }
    
    public void addGerp () {
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation/";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fList = Arrays.asList(fs);
    }
    
    public void addDAF () {
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation/";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String header = null;
            List<String> recordList = new ArrayList();
            String tem = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                header = br.readLine();
                List<String> l = PStringUtils.fastSplit(header);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tDAF\tDAF_ABD\t").append(l.get(9).replaceFirst("AAF", "DAF"));
                header = sb.toString();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextWriter(f.getAbsolutePath());
                bw.write(header);
                bw.newLine();
                float daf = -1;
                float dafABD = -1;
                float dafOther = -1;
                String subMajor = null;
                String subMinor = null;
                double subMaf = -1;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t");
                    l = PStringUtils.fastSplit(recordList.get(i));
                    if (l.get(5).equals(l.get(14))) {
                        sb.append((float)Double.parseDouble(l.get(7))).append("\t");
                    }
                    else if (l.get(6).equals(l.get(14))) {
                        sb.append((float)(1- Double.parseDouble(l.get(7)))).append("\t");
                    }
                    else sb.append("NA\t");
                    if (Double.parseDouble(l.get(8)) < 0.5) {
                        subMajor = l.get(3);
                        subMinor = l.get(4);
                        subMaf = Double.parseDouble(l.get(8));
                    }
                    else {
                        subMajor = l.get(4);
                        subMinor = l.get(3);
                        subMaf = 1 - Double.parseDouble(l.get(8));
                    }
                    if (l.get(14).equals(subMajor)) {
                        sb.append((float)subMaf).append("\t");
                    }
                    else if (l.get(14).equals(subMinor)) {
                        sb.append((float)(1-subMaf)).append("\t");
                    }
                    else sb.append("NA\t");
                    
                    tem = recordList.get(i);
                    if (Double.parseDouble(l.get(9)) < 0.5) {
                        subMajor = l.get(3);
                        subMinor = l.get(4);
                        subMaf = Double.parseDouble(l.get(9));
                    }
                    else {
                        subMajor = l.get(4);
                        subMinor = l.get(3);
                        subMaf = 1 - Double.parseDouble(l.get(9));
                    }
                    if (l.get(14).equals(subMajor)) {
                        sb.append((float)subMaf);
                    }
                    else if (l.get(14).equals(subMinor)) {
                        sb.append((float)(1-subMaf));
                    }
                    else sb.append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                System.out.println(tem);
                e.printStackTrace();
            }
        });
    }
    
    public void addAncestral () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/002_ancestral/byChr";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation/";
        File[] fs = new File (inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String annoFileS = f.getName().split("_")[0]+"_SNP_anno.txt";
            annoFileS = new File(dirS, annoFileS).getAbsolutePath();
            String header = null;
            List<String> recordList = new ArrayList();
            TIntArrayList posList = new TIntArrayList();
            String[] ancestral = null;
            try {
                BufferedReader br = IOUtils.getTextReader(annoFileS);
                header = br.readLine();
                String temp = null;
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                    l = PStringUtils.fastSplit(temp);
                    posList.add(Integer.parseInt(l.get(2)));
                }
                ancestral = new String[posList.size()];
                for (int i = 0; i < ancestral.length; i++) ancestral[i] = "NA";
                br.close();
                br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    int index = posList.binarySearch(pos);
                    if (index<0) continue;
                    ancestral[index] = l.get(2);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextWriter(annoFileS);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tAncestral");
                bw.write(sb.toString());
                bw.newLine();
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t").append(ancestral[i]);
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
    
    public void addSift2 () {
        String siftDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/001_sift/output";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation";
        File[] fs = new File(siftDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".xls");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            //f = new File("/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/sift/output/chr001.subgenome.maf0.01byPop.SNP_SIFTannotations.xls");
            String dbFileS = f.getName().split("\\.")[0]+"_SNP_anno.txt";
            dbFileS = new File (dirS, dbFileS).getAbsolutePath();
            RowTable<String> t = new RowTable (f.getAbsolutePath());
            SIFTRecord[] records = new SIFTRecord[t.getRowNumber()];
            for (int i = 0; i < records.length; i++) {
                SIFTRecord s = new SIFTRecord(Integer.parseInt(t.getCell(i, 1)), t.getCell(i, 3), t.getCell(i, 4), t.getCell(i, 7), t.getCell(i, 8), t.getCell(i, 12)); 
                records[i] = s;
            }
            Arrays.sort(records);
            try {
                List<String> dbList = new ArrayList();
                String temp = null;
                BufferedReader br = IOUtils.getTextReader(dbFileS);
                String header = br.readLine();
                while ((temp = br.readLine()) != null) {
                    dbList.add(temp);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextWriter(dbFileS);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tRegion\tVariant_type\tSIFT_score");
                bw.write(sb.toString());
                bw.newLine();
                List<String> l = null;
                for (int i = 0; i < dbList.size(); i++) {
                    l = PStringUtils.fastSplit(dbList.get(i));
                    SIFTRecord query = new SIFTRecord(Integer.parseInt(l.get(2)), l.get(4), l.get(10));
                    int index = Arrays.binarySearch(records, query);
                    if (index < 0) continue;
                    sb.setLength(0);
                    sb.append(dbList.get(i)).append("\t");
                    sb.append(records[index].region).append("\t");
                    sb.append(records[index].type).append("\t");
                    sb.append(records[index].value);
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
    
    class SIFTRecord implements Comparable<SIFTRecord> {

        public int pos;
        public String alt;
        public String transcript;
        public String region;
        public String type;
        public String value;

        
        public SIFTRecord(int pos, String alt, String transcript) {
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
        }
        
        public SIFTRecord(int pos, String alt, String transcript, String region, String type, String value) {
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
            this.region = region;
            this.type = type;
            this.value = value;
        }

        @Override
        public int compareTo(SIFTRecord o) { //
            if (this.pos < o.pos) { 
                return -1;
                } else if (this.pos == o.pos) {
                    int index = this.alt.compareTo(o.alt);
                    if (index < 0) {
                        return -1;
                    }
                    else if (index > 0) {
                        return 1;
                    }
                    else return transcript.compareTo(o.transcript);
                } else {
                    return 1;
            }
        }
}
    /**
     * @deprecated 
     */
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
            String outfileS = new File (outDirS, f.getName().replaceFirst(".vcf.gz", "_genicSNP.txt.gz")).getAbsolutePath();
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
