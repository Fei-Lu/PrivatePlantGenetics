/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import format.table.ColumnTable;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.Histogram;
import utils.IOUtils;
import utils.PArrayUtils;
import utils.PStringUtils;
import utils.wheat.RefV1Utils;

/**
 *
 * @author feilu
 */
public class AnnoGerp {
    
    public AnnoGerp() {
        //this.splitByChrID();
        //this.gerpStats();
        //this.gerpDistribution();
        //this.mkGerpDistributionGraph();
    }

    public void mkGerpDistributionGraph () {
        String infileS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/gerpSample.txt";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/gerpSample.pdf";
        ColumnTable<String> t = new ColumnTable<>(infileS);
        double[] values = t.getColumnAsDoubleArray(0);
        values = PArrayUtils.getNonredundantRandomSubset(values, 20000);
        Histogram h = new Histogram (values);
        h.setTitle("Gerp value distribution of wheat by aligning 29 species");
        h.setXLab("Gerp value");
        h.setYLab("Proportion");
        h.saveGraph(outfileS);
    }

    public void gerpDistribution () {
        String gerpDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/byChr";
        List<File> fList = IOUtils.getFileListInDirEndsWith(gerpDirS, ".gz");
        String outfileS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/gerpSample.txt";
        TIntArrayList[] siteLists = new TIntArrayList[fList.size()];
        TFloatArrayList[] gerpLists = new TFloatArrayList[fList.size()];
        double thresh = 0.00001;
        int step = 20000;
        fList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().split("_")[0].replaceFirst("chr", ""))-1;
            siteLists[chrIndex] = new TIntArrayList();
            gerpLists[chrIndex] = new TFloatArrayList();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<String> l = null;
                float value = 0;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    //if (Math.random() > thresh) continue;
                    cnt++;
                    if (cnt%step != 0) continue;
                    l = PStringUtils.fastSplit(temp);
                    value = Float.parseFloat(l.get(2));
                    siteLists[chrIndex].add(Integer.parseInt(l.get(1)));
                    gerpLists[chrIndex].add(value);
                }
                br.close();
                System.out.println(f.getAbsolutePath());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            String header = "Gerp";
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            StringBuilder sb = new StringBuilder ();
            for (int i = 0; i < siteLists.length; i++) {
                for (int j = 0; j < siteLists[i].size(); j++) {
                    sb.setLength(0);
                    //sb.append(i+1).append("\t").append(siteLists[i].get(j)).append("\t").append(gerpLists[i].get(j));
                    sb.append(gerpLists[i].get(j));
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
    }

    public void gerpStats () {
        String gerpDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/byChr";
        List<File> fList = IOUtils.getFileListInDirEndsWith(gerpDirS, ".gz");
        String gerpStats = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/gerpStats.txt";
        String header = "ChrID\tSites\tSites(>0)\tSites(>1)\tAlignRatio\tAlignRatio(>0)\tAlignRatio(>1)";
        int[] siteCounts = new int[fList.size()];
        int[] siteL0Counts = new int[fList.size()];
        int[] siteL1Counts = new int[fList.size()];
        fList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().split("_")[0].replaceFirst("chr", ""))-1;
            int siteCount = 0;
            int siteL0Count = 0;
            int siteL1Count = 0;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<String> l = null;
                double value = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    siteCount++;
                    value = Double.parseDouble(l.get(2));
                    if (value > 0) siteL0Count++;
                    if (value > 1) siteL1Count++;
                }
                br.close();
                siteCounts[chrIndex] = siteCount;
                siteL0Counts[chrIndex] = siteL0Count;
                siteL1Counts[chrIndex] = siteL1Count;
                System.out.println(f.getAbsolutePath());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IOUtils.getTextWriter(gerpStats);
            bw.write(header);
            bw.newLine();
            StringBuilder sb = new StringBuilder ();
            double value = 0;
            for (int i = 0; i < siteCounts.length; i++) {
                sb.setLength(0);
                sb.append(i+1).append("\t").append(siteCounts[i]).append("\t").append(siteL0Counts[i]).append("\t").append(siteL1Counts[i]).append("\t");
                value = (double)siteCounts[i]/RefV1Utils.getChrIDLength(i+1);
                sb.append((float)value).append("\t");
                value = (double)siteL0Counts[i]/RefV1Utils.getChrIDLength(i+1);
                sb.append((float)value).append("\t");
                value = (double)siteL1Counts[i]/RefV1Utils.getChrIDLength(i+1);
                sb.append((float)value);
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

    public void splitByChrID () {
        String wigDirS = "/Volumes/Fei_HDD_Mac/Gerp/asAlle";
        String gerpDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/byChr";
        String header = "Chr\tPos\tGerp";
        File[] fs = new File (wigDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".wig");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int[] chrID = new int[2];
            BufferedWriter[] bws = new BufferedWriter[2]; 
            String chromosome = f.getName().replaceFirst(".gerp.wig", "");
            chrID[0] = RefV1Utils.getChrID(chromosome, 1);
            chrID[1] = chrID[0]+1;
            try {
                for (int i = 0; i < chrID.length; i++) {
                    String outfileS = new File(gerpDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID[i])+"_gerp.txt.gz").getAbsolutePath();
                    bws[i] = IOUtils.getTextGzipWriter(outfileS);
                    bws[i].write(header);
                    bws[i].newLine();
                }
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = null;
                List<String> l = null;
                int chromPos = -1;
                int pos = -1;
                int index = -1;
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("f")) {
                        l = PStringUtils.fastSplitOnWhitespace(temp);
                        chromPos = Integer.parseInt(l.get(2).replaceFirst("start=", ""));
                    }
                    else {
                        sb.setLength(0);
                        pos = RefV1Utils.getPosOnChrID(chromosome, chromPos);
                        index = Arrays.binarySearch(chrID, RefV1Utils.getChrID(chromosome, chromPos));
                        sb.append(chrID[index]).append("\t").append(pos).append("\t").append(temp);
                        bws[index].write(sb.toString());
                        bws[index].newLine();
                        chromPos++;
                    }
                }
                for (int i = 0; i < chrID.length; i++) {
                    bws[i].flush();;
                    bws[i].close();
                }
                br.close();
                System.out.println(f.getAbsolutePath());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
}
