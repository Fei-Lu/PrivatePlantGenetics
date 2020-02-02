/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deprecated.analysis.maizeGeneticLoad;

import com.koloboke.collect.map.hash.HashIntFloatMap;
import com.koloboke.collect.map.hash.HashIntFloatMaps;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.format.genomeAnnotation.GeneFeature;
import pgl.format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import pgl.utils.IOFileFormat;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

/**
 *
 * @author feilu
 */
class VariantsAnnotation {
    
    public VariantsAnnotation () {
        //this.mkHeaderOfHap3();
        
        //this.convertAGPV3Annotation_deprecated();
        
        //this.convertAGPV3Anno();
        //this.sortAGPV4Anno();
        //this.mergeToRovertHmp321AGPV4();
        //this.addGerp();
    }
    
    public void addGerp () {
        String inputDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/hmp321_annotation/";
        String gerpDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/gerp/agpV4";
        String outputDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/annoDB";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        
        for (int i = 0; i < fs.length; i++) {
            //hmp321Info_chr010_AGPv4_Anno.txt
            int currentChr = Integer.valueOf(fs[i].getName().replaceFirst("hmp321Info_chr", "").replaceFirst("_AGPv4_Anno.txt", ""));
            String gerpFileS = new File (gerpDirS, "chr"+PStringUtils.getNDigitNumber(3, currentChr)+"_Gerp.txt").getAbsolutePath();
            String outputFileS = new File (outputDirS, "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, currentChr)+"_AGPv4_AnnoDB.txt").getAbsolutePath();
            String temp = null;
            try {
                BufferedReader br = IOUtils.getTextReader(gerpFileS);
                String header = br.readLine();
                List<String> l = PStringUtils.fastSplit(header);
                int cnt = 0;
                HashIntFloatMap posLengthMap = HashIntFloatMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
                HashIntFloatMap posValueMap = HashIntFloatMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+gerpFileS);
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    posLengthMap.put(pos, Float.parseFloat(l.get(2)));
                    posValueMap.put(pos, Float.parseFloat(l.get(3)));
                }
                br.close();
                //Collections.sort(rl);
                br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
                header  = br.readLine()+"\tGerpTreeLength\tGerpValue";
                bw.write(header);
                bw.newLine();
                int index = -1;
                StringBuilder sb = null;
                float len = -1;
                float value = -1;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    len = posLengthMap.get(pos);
                    value = posValueMap.get(pos);
                    sb = new StringBuilder(temp);
                    if (len == -1) {
                        sb.append("\tNA\tNA");
                    }
                    else {
                        sb.append("\t").append(len).append("\t").append(value);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        
    }
    
    public void mergeToRovertHmp321AGPV4 () {
        String inputDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/hmp321_header";
        String inputAnnoFileS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/source/003_hmp321Info_AGPv4";
        String outputDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/hmp321_annotation/";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        for (int i = 0; i < fs.length; i++) {
            //hmp321_agpv4_chr1.vcf.header.txt
            int currentChr = Integer.valueOf(fs[i].getName().replaceFirst("hmp321_agpv4_chr", "").replaceFirst(".vcf.header.txt", ""));
            //hmp321Info_chr001_AGPv4.txt
            String inAnnoS = new File (inputAnnoFileS, "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, currentChr)+"_AGPv4.txt").getAbsolutePath();
            String outAnnoS = new File (outputDirS, "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, currentChr)+"_AGPv4_Anno.txt").getAbsolutePath();
            String temp = null;
            try {
                List<Record> rl = new ArrayList<>();
                BufferedReader br = IOUtils.getTextReader(inAnnoS);
                String header = br.readLine();
                
                List<String> l = PStringUtils.fastSplit(header);
                int columnNum = l.size();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    if (currentChr != Integer.parseInt(l.get(0))) continue;
                    Record r = new Record (Integer.parseInt(l.get(1)), temp);
                    rl.add(r);
                }
                br.close();
                Collections.sort(rl);
                br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                temp = br.readLine();
                BufferedWriter bw = IOUtils.getTextWriter(outAnnoS);
                bw.write(header);
                bw.newLine();
                int index = -1;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    Record query = new Record (Integer.parseInt(l.get(1)), "");
                    index = Collections.binarySearch(rl, query);
                    if (index < 0) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3));
                        for (int j = 0; j < columnNum-4; j++) {
                            sb.append("\t").append("NA");
                        }
                        bw.write(sb.toString());
                    }
                    else {
                        bw.write(rl.get(index).r);
                    }
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                System.out.println(temp);
                e.printStackTrace();
            }
        }
    }
    
    public void sortAGPV4Anno () {
        String infileDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/source/002_hmp321Info_AGPv4";
        String outputDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/source/003_hmp321Info_AGPv4";
        File[] fs = new File (infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".tem");
        int chrNum = fs.length;
        
        for (int i = 0; i < chrNum; i++) {
            int currentChr = i+1;
            List<Record> rl = new ArrayList<>();
            String header = null;
            String outfileS = new File (outputDirS, "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, currentChr)+"_AGPv4.txt").getAbsolutePath();
            for (int j = 0; j < fs.length; j++) {
                try {
                    BufferedReader br = IOUtils.getTextReader(fs[j].getAbsolutePath());
                    header = br.readLine();
                    String temp = null;
                    List<String> l = null;
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp.substring(0, 15));
                        if (currentChr != Integer.parseInt(l.get(0))) continue;
                        Record r = new Record (Integer.parseInt(l.get(1)), temp);
                        rl.add(r);
                    }
                    br.close();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            Collections.sort(rl);
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(header);
                bw.newLine();
                for (int j = 0; j < rl.size(); j++) {
                    bw.write(rl.get(j).r);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        
    }
    
    class Record implements Comparable<Record> {
        int pos;
        public String r;
        public Record (int pos, String r) {
            this.pos = pos;
            this.r = r;
        }
        @Override
        public int compareTo(Record o) {
            return this.pos-o.pos;
        }
    }
    
    public void convertAGPV3Anno () {
        String infileDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/source/001_hmp321Info/";
        String outfileDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/source/002_hmp321Info_AGPv4/";
        String tempDirS = "/Users/feilu/Documents/database/maize/crossMap/test";
        File[] fs = new File (infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
            String outfileS = new File (outfileDirS, f.getName().replaceFirst(".gz", ".tem")).getAbsolutePath();
            String tempFileS = new File(tempDirS, f.getName().replaceFirst(".gz", ".bed")).getAbsolutePath();
            RowTable<String> t = new RowTable (f.getAbsolutePath());
            int[] chr = new int[t.getRowNumber()];
            int[] pos = new int[t.getRowNumber()];
            for (int i = 0; i < chr.length; i++) {
                chr[i] = Integer.parseInt(t.getCell(i, 0));
                pos[i] = Integer.parseInt(t.getCell(i, 1));
            }
            CrossMapUtils cm = new CrossMapUtils(chr, pos, tempFileS);
            cm.setMaizeV3ToV4();
            cm.convert();
            List<int[]> l = cm.getConvertedCoordinate();
            cm.deleteBedFiles();
            chr = l.get(0);
            pos = l.get(1);
            for (int i = 0; i < t.getRowNumber(); i++) {
                t.setCell(i, 0, String.valueOf(chr[i]));
                t.setCell(i, 1, String.valueOf(pos[i]));
            }
            t.writeTextTable(outfileS, IOFileFormat.Text);
        });
    }
    
    /**
     * @deprecated 
     */
    public void convertAGPV3Annotation_deprecated () {
        String v3AnnotationDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/source/001_hmp321Info";
        String v3V4MapDirS = "/Users/feilu/Documents/database/maize/crossMap/V3V4Map";
        String outputDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/source/002_hmp321Info_AGPv4";
        File[] fs = new File (v3V4MapDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".map");
        List<File> fList = Arrays.asList(fs);
        HashIntIntMap[] v3v4Maps = new HashIntIntMap[fList.size()];
        fList.parallelStream().forEach(f -> {
            int index = Integer.valueOf(f.getName().replaceFirst(".map", "").replaceFirst("chr", ""))-1;
            v3v4Maps[index] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    l = PStringUtils.fastSplit(temp);
                    if (l.get(2).startsWith("N")) continue;
                    v3v4Maps[index].put(Integer.parseInt(l.get(1)), Integer.parseInt(l.get(3)));
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        fs = new File(v3AnnotationDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int index = Integer.parseInt(f.getName().replaceFirst("hmp321Info_chr", "").replaceFirst(".txt.gz", ""))-1;
            String outfileS = new File (outputDirS, "AGPv4_"+f.getName()).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                String header = temp;
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                List<String> l = null;
                StringBuilder sb = null;
                bw.write(header);
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    sb = new StringBuilder();
                    l = PStringUtils.fastSplit(temp);
                    int query = Integer.parseInt(l.get(1));
                    int pos = v3v4Maps[index].get(query);
                    if (pos == -1) {
                        sb.append("NA\tNA");
                    }
                    else {
                        sb.append(l.get(0)).append("\t").append(pos);
                    }
                    for (int i = 2; i < l.size(); i++) {
                        sb.append("\t").append(l.get(i));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void mkHeaderOfHap3 () {
        String infileDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/hmp321_agp4";
        String outfileDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/hmp321_header";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String oufileS = new File(outfileDirS, f.getName().replaceFirst(".gz", ".header.txt")).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(oufileS);
                bw.write("Chr\tPos\tRef\tAlt\tMaf");
                bw.newLine();
                String temp = null;
                List<String> l = null;
                StringBuilder sb = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    cnt++;
                    if (cnt%100000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    sb = new StringBuilder();
                    temp = temp.substring(0, 160);
                    l = PStringUtils.fastSplit(temp);
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    l = PStringUtils.fastSplit(l.get(7), "MAF=");
                    l = PStringUtils.fastSplit(l.get(1), ";");
                    sb.append(l.get(0));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
}
