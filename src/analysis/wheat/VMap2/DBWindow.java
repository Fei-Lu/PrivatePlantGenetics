package analysis.wheat.VMap2;

import pgl.infra.genomeAnnotation.GeneFeature;
import pgl.infra.range.Range;
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
import java.util.List;

public class DBWindow {
    int windowSize = 2000000;
    int windowStep = 1000000;

    public DBWindow () {
        //this.initializeWithDelAndSyn();
        //this.addCDSLengthInWindow();
    }

    public void addCDSLengthInWindow () {
        String geneFeatureFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String dbFileS = "/Users/feilu/Documents/analysisH/vmap2/005_genomeWindowDB/windowDB.txt";
        String hcGeneFileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        RowTable<String> gt = new RowTable<>(hcGeneFileS);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();
        List<String> chromosomeList = RefV1Utils.getChromosomeList();
        TIntArrayList cdsWindowList = new TIntArrayList();
        for (int i = 0; i < chromosomeList.size(); i++) {
            int chrlength = RefV1Utils.getChromosomeLength(chromosomeList.get(i));
            SimpleWindow sw = new SimpleWindow(chrlength, windowSize, windowStep);
            int chrID = RefV1Utils.getChrID(chromosomeList.get(i), 1);
            int geneIndex = -1;
            int tranIndex = -1;
            int cdsStart = -1;
            int cdsEnd = -1;
            for (int j = 0; j < gt.getRowNumber(); j++) {
                int currentChrID = Integer.parseInt(gt.getCell(j, 2));
                if (currentChrID < chrID) continue;
                else if (currentChrID > (chrID +1)) break;
                geneIndex = gf.getGeneIndex(gt.getCell(j, 0));
                for (int k = 0; k < gf.getTranscriptNumber(geneIndex); k++) {
                    if (!gt.getCell(j,1).equals(gf.getTranscriptName(geneIndex, k)))continue;
                    tranIndex = k;
                    break;
                }
                List<Range> cdsList = gf.getCDSList(geneIndex, tranIndex);
                for (int k = 0; k < cdsList.size(); k++) {
                    cdsStart = cdsList.get(k).start;
                    cdsEnd = cdsList.get(k).end;
                    if (currentChrID == chrID + 1) {
                        cdsStart = RefV1Utils.getPosOnChromosome(currentChrID, cdsStart);
                        cdsEnd = RefV1Utils.getPosOnChromosome(currentChrID, cdsEnd);
                    }
                    sw.addPositionCountFromRange(cdsStart, cdsEnd);
                }
            }
            cdsWindowList.add(sw.getWindowValuesInt());
        }
        Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(dbFileS);
        String header = two.getFirstElement();
        List<String> recordList = two.getSecondElement();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(dbFileS);
            header = header + "\tCDSLength\tDelFrequency\tSynFrequency";
            bw.write(header);
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            List<String> l = new ArrayList<>();
            int cdsLength = -1;
            for (int i = 0; i < recordList.size(); i++) {
                sb.setLength(0);
                l = PStringUtils.fastSplit(recordList.get(i));
                cdsLength = cdsWindowList.get(i);
                sb.append(recordList.get(i)).append("\t").append(cdsLength).append("\t");
                if (cdsLength != 0) {
                    sb.append((float)((double)Integer.parseInt(l.get(3))/cdsLength)).append("\t");
                    sb.append((float)((double)Integer.parseInt(l.get(4))/cdsLength));
                }
                else {
                    sb.append(Float.NaN).append("\t").append(Float.NaN);
                }
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

    public void initializeWithDelAndSyn () {
        String delInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/del";
        String synInfoDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousBiology/001_snp/syn";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap2/005_genomeWindowDB/windowDB.txt";
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
                    String delFileS = new File(delInfoDirS, "chr"+ PStringUtils.getNDigitNumber(3, chrID)+"_SNP_anno.txt.gz").getAbsolutePath();
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
}
