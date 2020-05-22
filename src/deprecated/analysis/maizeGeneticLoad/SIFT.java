/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deprecated.analysis.maizeGeneticLoad;

import pgl.infra.anno.gene.GeneFeature;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author feilu
 */
class SIFT {
    
    public SIFT () {
        //this.testData();
        this.addSIFT();
    }
    
    public void addSIFT () {
        String infileS = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.pgf";
        String siftDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/sift/source/";
        String annotationDBDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/annoDB";
        String outputDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/annoDB_new";
        new File (outputDirS).mkdir();
        GeneFeature gf = new GeneFeature(infileS);
        int geneNum = gf.getGeneNumber();
        String[] trans = new String[geneNum];
        for (int i = 0; i < geneNum; i++) {
            int index = gf.getLongestTranscriptIndex(i);
            trans[i] = gf.getTranscriptName(i, index);
        }
        Arrays.sort(trans);
        File[] fs = new File(siftDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "xls");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().split("_")[2].replaceFirst("chr", ""))-1;
//            String infileSS = "/Users/feilu/Documents/analysisL/production/maizeLoad/sift/source/hmp321_agpv4_chr10_SIFTannotations.xls";
//            f = new File(infileSS);
//            chrIndex = 9;
            //hmp321Info_chr001_AGPv4_AnnoDB.txt
            String dbFileS = new File (annotationDBDirS, "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, chrIndex+1)+"_AGPv4_AnnoDB.txt").getAbsolutePath();
            String outfileS = new File(outputDirS,  "hmp321Info_chr"+PStringUtils.getNDigitNumber(3, chrIndex+1)+"_AGPv4_AnnoDB.txt").getAbsolutePath();
            String header = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<SIFTRecord> sList = new ArrayList<>();
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String alt = l.get(3);
                    String transcript = l.get(4).replaceFirst("transcript.", "");
                    String type = l.get(8);
                    String value = l.get(12);
                    SIFTRecord s = new SIFTRecord( pos, alt,  transcript,  type, value);
                    sList.add(s);
                }
                br.close();
                Collections.sort(sList);
                br = IOUtils.getTextReader(dbFileS);
                header = br.readLine()+"\tVariant_type\tSIFT_score\tTranscript";
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(header);
                bw.newLine();
                int index = -1;
                int startIndex = -1;
                int endIndex = -1;
                StringBuilder sb = null;
                while ((temp = br.readLine()) != null) {
                    sb = new StringBuilder();
                    sb.append(temp).append("\t");
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String alt = l.get(3).split(",")[0];
                    SIFTRecord query = new SIFTRecord (pos, alt, "", "", "");
                    index = Collections.binarySearch(sList, query);
                    if (index < 0) {
                        
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    startIndex = index;
                    endIndex = index;
                    while ((startIndex-1)>-1 && sList.get(startIndex-1).isSimilar(pos, alt)) {
                        startIndex--;
                    }
                    while ((endIndex+1) < sList.size() && sList.get(endIndex+1).isSimilar(pos, alt)) {
                        endIndex++;
                    }
                    boolean status = false;
                    for (int i = startIndex; i < endIndex+1; i++) {
                        if (Arrays.binarySearch(trans, sList.get(i).transcript) >= 0 && sList.get(i).type.equals("NONSYNONYMOUS") || sList.get(i).type.equals("SYNONYMOUS")) {
                            if (Arrays.binarySearch(trans, sList.get(i).transcript) < 0) {
                                sb.append("NA\tNA\tNA");
                            }
                            else {
                                sb.append(sList.get(i).type).append("\t").append(sList.get(i).value).append("\t").append(sList.get(i).transcript);
                            }
                            
                            bw.write(sb.toString());
                            bw.newLine();
                            status = true;
                            break;
                        }
                    }
                    if (status ==false) {
                        sb.append("NA\tNA\tNA");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
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
    
    class SIFTRecord implements Comparable<SIFTRecord>{
        public int pos;
        public String alt;
        public String transcript;
        public String type;
        public String value;
        
        public SIFTRecord (int pos, String alt, String transcript, String type, String value) {            
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
            this.type = type;
            this.value = value;
        }
        
        public boolean isSimilar (int pos, String alt) {
            if (pos == this.pos && alt.equals(this.alt)) return true;
            return false;
        }
        
        @Override
        public int compareTo(SIFTRecord o) {
            if (this.pos < o.pos) {
                return -1;
            }
            else if (this.pos == o.pos) {
                return this.alt.compareTo(o.alt);
            }
            else {
                return 1;
            }
            
        }
    }
    
    public boolean ifMatch (List<String> l, int pos, String alt, String[] trans) {
        int posS = Integer.parseInt(l.get(1));
        if (posS != pos) return false;
        if (!alt.equals(l.get(3))) return false;
        String type = l.get(8);
        if (!(type.startsWith("NONSYNONYMOUS") || type.startsWith("SYNONYMOUS"))) return false;
        String query = l.get(4).replaceFirst("transcript.", "");
        if (Arrays.binarySearch(trans, query) < 0) return false;
        return true;
    }
    
    public void testData () {
        String infileS = "/Users/feilu/Documents/analysisL/production/maizeLoad/sift/source/hmp321_agpv4_chr10_SIFTannotations.xls";
        Set<String> typeSet = new HashSet<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                typeSet.add(tem[8]);
            }
            br.close();
            String[] types = typeSet.toArray(new String[typeSet.size()]);
            Arrays.sort(types);
            Set<String>[] subsets = new HashSet[types.length];
            for (int i = 0; i < subsets.length; i++) subsets[i] = new HashSet<>();
            br = IOUtils.getTextReader(infileS);
            temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                String query = tem[8];
                int index = Arrays.binarySearch(types, query);
                subsets[index].add(tem[16]);
            }
            br.close();
            for (int i = 0; i < types.length; i++) {
                System.out.println(types[i]);
                System.out.println(subsets[i]+"\n");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
}
