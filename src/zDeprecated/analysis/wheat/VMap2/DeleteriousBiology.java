/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zDeprecated.analysis.wheat.VMap2;

import pgl.infra.table.ColumnTable;
import pgl.infra.table.RowTable;
import pgl.infra.window.SimpleWindow;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;

/**
 *
 * @author feilu
 */
public class DeleteriousBiology {
    
    public DeleteriousBiology () {
        //this.countDeleteriousAndSyn();
        //this.identifyDeleteriousAndSyn();
        //this.mkVCFofDeleteriousAndSyn();
        //this.chromosomeDistribution();
    }

    /**
     * @deprecated
     */
    public void chromosomeDistribution () {
        int windowSize = 20000000;
        int windowStep = 5000000;
        String delInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/del";
        String synInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/syn";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/003_chrDis/delSynOnChr.txt";
        String header = "Chromosome\tWindowStart\tWindowEnd\tDelCount\tSynCount\tDelSynRatio";
        List<String> chromList = RefV1Utils.getChromosomeList();
        StringBuilder sb = new StringBuilder();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            List<String> l = new ArrayList<>();
            for (int i = 0; i < chromList.size(); i++) {
                int chrLength = RefV1Utils.getChromosomeLength(chromList.get(i));
                TIntArrayList delList = new TIntArrayList();
                TIntArrayList synList = new TIntArrayList();
                int chrID = RefV1Utils.getChrID(chromList.get(i), 1);
                for (int j = 0; j < 2; j++) {
                    chrID+=j;
                    String delFileS = new File (delInfoDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_SNP_anno.txt.gz").getAbsolutePath();
                    String synFileS = new File (synInfoDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_SNP_anno.txt.gz").getAbsolutePath();
                    BufferedReader br = IOUtils.getTextGzipReader(delFileS);
                    String temp = br.readLine();
                    int pos = -1;
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp.substring(0, 50));
                        pos = Integer.parseInt(l.get(2));
                        if (j == 1) pos = RefV1Utils.getPosOnChromosome(chrID, pos);
                        delList.add(pos);
                    }
                    br.close();
                    br = IOUtils.getTextGzipReader(synFileS);
                    temp = br.readLine();
                    pos = -1;
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp.substring(0, 50));
                        pos = Integer.parseInt(l.get(2));
                        if (j == 1) pos = RefV1Utils.getPosOnChromosome(chrID, pos);
                        synList.add(pos);
                    }
                }
                SimpleWindow sw = new SimpleWindow(chrLength, windowSize, windowStep);
                sw.addPositionCount(delList.toArray());
                int[] delWindowCount = sw.getWindowValuesInt();
                sw.clearWindowValues();
                sw.addPositionCount(synList.toArray());
                int[] synWindowCount = sw.getWindowValuesInt();
                int[] windowStarts = sw.getWindowStarts();
                int[] windowEnds = sw.getWindowEnds();
                for (int j = 0; j < windowStarts.length; j++) {
                    sb.setLength(0);
                    sb.append(chromList.get(i)).append("\t").append(windowStarts[j]).append("\t").append(windowEnds[j]).append("\t");
                    sb.append(delWindowCount[j]).append("\t").append(synWindowCount[j]).append("\t").append((float)((double)delWindowCount[j]/synWindowCount[j]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(chromList.get(i));
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void mkVCFofDeleteriousAndSyn () {
        String vmapDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/002_exonSNPVCF/";
        String delInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/del";
        String synInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/syn";
        String delVcfDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/002_vcf/del";
        String synVcfDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/002_vcf/syn";
        List<File> fList = IOUtils.getFileListInDirEndsWith(vmapDirS, ".gz");
        fList.parallelStream().forEach(f -> {
            String chrS = f.getName().split("_")[0];
            String delInfoFileS = new File (delInfoDirS, chrS+"_SNP_anno.txt.gz").getAbsolutePath();
            String synInfoFileS = new File (synInfoDirS, chrS+"_SNP_anno.txt.gz").getAbsolutePath();
            String delVcfFileS = new File (delVcfDirS, chrS+"_del_vmap2.1.vcf.gz").getAbsolutePath();
            String synVcfFileS = new File (synVcfDirS, chrS+"_syn_vmap2.1.vcf.gz").getAbsolutePath();
            try {
                ColumnTable<String> t = new ColumnTable<>(delInfoFileS);
                int[] delsites = t.getColumnAsIntArray(2);
                Arrays.sort(delsites);
                t = new ColumnTable<>(synInfoFileS);
                int[] synsites = t.getColumnAsIntArray(2);
                Arrays.sort(synsites);
                BufferedWriter bwd = IOUtils.getTextGzipWriter(delVcfFileS);
                BufferedWriter bws = IOUtils.getTextGzipWriter(synVcfFileS);
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()).startsWith("##")) {
                    bwd.write(temp); bwd.newLine();
                    bws.write(temp); bws.newLine();
                }
                bwd.write(temp); bwd.newLine();
                bws.write(temp); bws.newLine();
                List<String> l = new ArrayList<>();
                int index = -1;
                int pos = -1;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp.substring(0, 60));
                    pos = Integer.parseInt(l.get(1));
                    if (Arrays.binarySearch(delsites, pos) >= 0) {
                        bwd.write(temp);bwd.newLine();
                    }
                    else if (Arrays.binarySearch(synsites, pos) >= 0) {
                        bws.write(temp);bws.newLine();
                    }
                }
                bwd.flush();bwd.close();
                bws.flush();bws.close();
                br.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void identifyDeleteriousAndSyn () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        String delDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/del";
        String synDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/syn";
        double siftThresh = 0.05;
        double gerpThresh = 1;
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        LongAdder counterDel = new LongAdder();
        LongAdder counterSyn = new LongAdder();
        fList.stream().forEach(f -> {
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String delFileS = new File (delDirS, f.getName()).getAbsolutePath();
            String synFileS = new File (synDirS, f.getName()).getAbsolutePath();
            String header = two.getFirstElement();
            List<String> recordList = two.getSecondElement();
            int typeIndex = -1;
            int siftIndex = -1;
            int gerpIndex = -1;
            int cnt = 0;
            int cntSyn = 0;
            String[] temp = header.split("\t");
            for (int i = 0; i < temp.length; i++) {
                if (temp[i].equals("Variant_type")) typeIndex = i;
                else if (temp[i].equals("Derived_SIFT")) siftIndex = i;
                else if (temp[i].equals("Gerp")) gerpIndex = i;
            }
            List<String> l = null;
            try {
                BufferedWriter bwd = IOUtils.getTextGzipWriter(delFileS);
                BufferedWriter bws = IOUtils.getTextGzipWriter(synFileS);
                bwd.write(header);
                bwd.newLine();
                bws.write(header);
                bws.newLine();
                for (int i = 0; i < recordList.size(); i++) {
                    l = PStringUtils.fastSplit(recordList.get(i));
                    if (l.get(typeIndex).startsWith("SYN")) {
                        cntSyn++;
                        bws.write(recordList.get(i));
                        bws.newLine();
                    }
                    if (!l.get(typeIndex).startsWith("NONS")) continue;
                    if (l.get(siftIndex).startsWith("N")) continue;
                    if (l.get(gerpIndex).startsWith("N")) continue;
                    if (Double.parseDouble(l.get(siftIndex)) < siftThresh && Double.parseDouble(l.get(gerpIndex)) > gerpThresh) {
                        cnt++;
                        bwd.write(recordList.get(i));
                        bwd.newLine();
                    }
                }
                bwd.flush();
                bwd.close();
                bws.flush();
                bws.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            StringBuilder sb = new StringBuilder();
            sb.append("Del:\t").append(cnt).append("\tSyn:\t").append(cntSyn).append("\t").append(f.getName());
            System.out.println(sb.toString());
            counterDel.add(cnt);
            counterSyn.add(cntSyn);
        });
        StringBuilder sb = new StringBuilder();
        sb.append("Del:\t").append(counterDel.intValue()).append("\t").append(counterSyn.intValue());
        System.out.println(sb.toString());
    }

    public void countDeleteriousAndSyn() {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        double siftThresh = 0.05;
        double gerpThresh = 1;
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        LongAdder counterDel = new LongAdder();
        LongAdder counterSyn = new LongAdder();
        fList.stream().forEach(f -> {
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            int typeIndex = t.getColumnIndex("Variant_type");
            int siftIndex = t.getColumnIndex("Derived_SIFT");
            int gerpIndex = t.getColumnIndex("Gerp");
            int cnt = 0;
            int cntSyn = 0;
            for (int i = 0; i < t.getRowNumber(); i++) {
                if (t.getCell(i, typeIndex).startsWith("SYN")) cntSyn++;
                if (!t.getCell(i, typeIndex).startsWith("NONS")) continue;
                if (t.getCell(i, siftIndex).startsWith("N")) {
                    continue;
                }
                if (t.getCell(i, gerpIndex).startsWith("N")) continue;
                if (Double.parseDouble(t.getCell(i, siftIndex)) < siftThresh && Double.parseDouble(t.getCell(i, gerpIndex)) > gerpThresh) {
                    cnt++;
                }
            }
            StringBuilder sb = new StringBuilder();
            sb.append("Del:\t").append(cnt).append("\tSyn:\t").append(cntSyn).append("\t").append(f.getAbsoluteFile());
            System.out.println(sb.toString());
            counterDel.add(cnt);
            counterSyn.add(cntSyn);
        });
        StringBuilder sb = new StringBuilder();
        sb.append("Del:\t").append(counterDel.intValue()).append("\t").append(counterSyn.intValue());
        System.out.println(sb.toString());
    }
}
