package analysis.wheat.VMap1;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.TIntHashSet;
import htsjdk.samtools.util.IOUtil;
import pgl.graph.r.DensityPlot;
import pgl.graph.r.Histogram;
import pgl.infra.range.Range;
import pgl.infra.range.Ranges;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

class Introgression {

    public Introgression () {
//        this.intervalSize();
//        this.findIndividualWithMaxFd();
//        this.maxFdHist();
    }


    public void maxFdHist () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/maxFd";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/maxFd_hist";
        List<File> dirList = IOUtils.getDirListInDir(inDirS);
        for (int i = 0; i < dirList.size(); i++) {
            String outfileS = new File (outDirS,dirList.get(i).getName()+".pdf").getAbsolutePath();
            List<File> fList = IOUtils.getFileListInDirEndsWith(dirList.get(i).getAbsolutePath(), ".gz");
            TDoubleArrayList fds = new TDoubleArrayList();
            for (int j = 0; j < fList.size(); j++) {
                try {
                    BufferedReader br = IOUtils.getTextGzipReader(fList.get(j).getAbsolutePath());
                    String temp = br.readLine();
                    List<String> l = new ArrayList<>();
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp);
                        fds.add(Double.parseDouble(l.get(3)));
                    }
                    br.close();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            System.out.println(outfileS);
            double[] sample = PArrayUtils.getNonredundantRandomSubset(fds.toArray(), 20000);
            Histogram h = new Histogram(sample);
            h.setTitle(dirList.get(i).getName());
            h.setXLab("fd");
            h.setYLab("Proportion");
            h.setBreakNumber(100);
            h.setXLim(0,1);
            h.saveGraph(outfileS);
        }
    }

    public void findIndividualWithMaxFd () {
        String inDirS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw";
        List<File> dirListA = IOUtils.getDirListInDirStartsWith(inDirS, "A");
        List<File> dirListB = IOUtils.getDirListInDirStartsWith(inDirS, "B");
        List<File> dirListD = IOUtils.getDirListInDirStartsWith(inDirS, "D");
        dirListA.addAll(dirListB);
        String sampleAFileS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/A025/A1.csv.gz";
        String sampleBFileS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/B001/AB1.csv.gz";
        String sampleDFileS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/D001/D1.csv.gz";
        String outfileDirSA = "/Users/feilu/Documents/analysisH/vmap1/fd/maxFd/A_maxFd";
        String outfileDirSB = "/Users/feilu/Documents/analysisH/vmap1/fd/maxFd/B_maxFd";
        String outfileDirSD = "/Users/feilu/Documents/analysisH/vmap1/fd/maxFd/D_maxFd";
        Ranges ranA = this.getRanges(sampleAFileS, "A");
        Ranges ranB = this.getRanges(sampleBFileS, "B");
        Ranges ranD = this.getRanges(sampleDFileS, "D");
        String[] taxaA = this.getTaxaNames(dirListA.get(0));
        String[] taxaB = this.getTaxaNames(dirListB.get(0));
        String[] taxaD = this.getTaxaNames(dirListD.get(0));

        this.outputMaxFd(dirListA, ranA, taxaA, outfileDirSA);
        this.outputMaxFd(dirListB, ranB, taxaB, outfileDirSB);
        this.outputMaxFd(dirListD, ranD, taxaD, outfileDirSD);
    }

    private void outputMaxFd (List<File> sourceDirList, Ranges rs, String[] targetTaxa, String outDirS) {
        File outDir = new File (outDirS);
        outDir.mkdir();
        rs.sortByStartPosition();
        for (int i = 0; i < targetTaxa.length; i++) {
            File[] inputFiles = new File[sourceDirList.size()];
            float[][] fds = new float[rs.getRangeNumber()][sourceDirList.size()];
            RowTable<String> t = null;
            for (int j = 0; j < sourceDirList.size(); j++) {
                if (sourceDirList.get(j).getName().startsWith("B")) {
                    inputFiles[j] = new File(sourceDirList.get(j), "AB"+targetTaxa[i]+".csv.gz").getAbsoluteFile();
                }
                else if (sourceDirList.get(j).getName().startsWith("A")) {
                    inputFiles[j] = new File(sourceDirList.get(j), "A"+targetTaxa[i]+".csv.gz").getAbsoluteFile();
                }
                else if (sourceDirList.get(j).getName().startsWith("D")) {
                    inputFiles[j] = new File(sourceDirList.get(j), "D"+targetTaxa[i]+".csv.gz").getAbsoluteFile();
                }
                t = new RowTable<>(inputFiles[j].getAbsolutePath(), ",");
                for (int k = 0; k < t.getRowNumber(); k++) {
                    Range r = new Range (Integer.parseInt(t.getCell(k, 0)), Integer.parseInt(t.getCell(k,1)), Integer.parseInt(t.getCell(k,2)));
                    int index = Collections.binarySearch(rs.getRangeList(), r);
                    if (index < 0) continue;
                    String dS = t.getCell(k,8);
                    String fdS = t.getCell(k, 9);
                    if (dS.startsWith("n") || Double.parseDouble(dS) < 0) {
                        fds[index][j] = 0;
                    }
                    else if (fdS.startsWith("n")) {
                        fds[index][j] = 0;
                    }
                    else {
                        float v = Float.parseFloat(fdS);
                        if (v < 0 || v > 1) fds[index][j] = 0;
                        else fds[index][j] = (float)v;
                    }
                }
            }
            String outfileS = new File (outDirS, targetTaxa[i]+"_maxFd.txt.gz").getAbsolutePath();
            int cnt = 0;
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("Chr\tStart\tEnd\tMaxFd\tTaxa");
                bw.newLine();
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < rs.getRangeNumber(); j++) {
                    sb.setLength(0);
                    sb.append(rs.getRangeChromosome(j)).append("\t").append(rs.getRangeStart(j)).append("\t").append(rs.getRangeEnd(j)).append("\t");
                    float max = 0;
                    for (int k = 0; k < fds[j].length; k++) {
                        if (fds[j][k] > max) max = fds[j][k];
                    }
                    sb.append(max);
                    for (int k = 0; k < fds[j].length; k++) {
                        cnt = k;
                        if (fds[j][k] == max) sb.append("\t").append(sourceDirList.get(k).getName());
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.out.println(cnt);
            }
        }
    }

    private String[] getTaxaNames (File dir) {
        List<File> fList = IOUtils.getFileListInDirEndsWith(dir.getAbsolutePath(), ".gz");
        String[] taxa = new String[fList.size()];
        for (int i = 0; i < taxa.length; i++) {
            taxa[i] = String.valueOf(i+1);
        }
        return taxa;
    }

    private Ranges getRanges (String sampleFileS, String subgenome) {
        RowTable<String> t = new RowTable<>(sampleFileS, ",");
        List<Range> rList = new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chr = Integer.parseInt(t.getCell(i,0));
            String type = RefV1Utils.getSubgenomeFromChrID(chr);
            if (type.equals(subgenome)) {
                Range r = new Range(chr, Integer.parseInt(t.getCell(i,1)), Integer.parseInt(t.getCell(i,2)));
                rList.add(r);
            }
        }
        return new Ranges(rList);
    }

    public void intervalSize() {
//        String file1 = "/Users/feilu/Desktop/untitled folder/A1.csv";
//        String file2 = "/Users/feilu/Desktop/untitled folder/A2.csv";
//        String file3 = "/Users/feilu/Desktop/untitled folder/AB1.csv";
//        String file4 = "/Users/feilu/Desktop/untitled folder/AB2.csv";
        String infileS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/B001/AB1.csv.gz";
        String infileSD = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/D001/D1.csv.gz";
        String intervalA = "/Users/feilu/Documents/analysisH/vmap1/fd/interval/A_interval_density.pdf";
        String intervalB = "/Users/feilu/Documents/analysisH/vmap1/fd/interval/B_interval_density.pdf";
        String intervalD = "/Users/feilu/Documents/analysisH/vmap1/fd/interval/D_interval_density.pdf";

        RowTable<String> t = new RowTable<>(infileS, ",");
        int[] chrs = t.getColumnAsIntArray(0);
        TIntHashSet cSet = new TIntHashSet(chrs);
        chrs = cSet.toArray();
        Arrays.sort(chrs);
        TDoubleArrayList AList = new TDoubleArrayList();
        TDoubleArrayList BList = new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String genomeType = RefV1Utils.getSubgenomeFromChrID(Integer.parseInt(t.getCell(i,0)));
            int dis = Integer.parseInt(t.getCell(i,2))-Integer.parseInt(t.getCell(i,1));
            if (genomeType.equals("A")) {
                AList.add(dis);
            }
            else if (genomeType.equals("B")) {
                BList.add(dis);
            }
        }

        t = new RowTable<>(infileSD, ",");
        chrs = t.getColumnAsIntArray(0);
        cSet = new TIntHashSet(chrs);
        chrs = cSet.toArray();
        Arrays.sort(chrs);
        TDoubleArrayList DList = new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String genomeType = RefV1Utils.getSubgenomeFromChrID(Integer.parseInt(t.getCell(i,0)));
            int dis = Integer.parseInt(t.getCell(i,2))-Integer.parseInt(t.getCell(i,1));
            if (genomeType.equals("D")) {
                DList.add(dis);
            }
        }


        DensityPlot d = new DensityPlot(AList.toArray());
        d.setSmoothN(5000);
        d.setTitle("Interval distribution of fd test in A subgenome");
        d.setXLab("Interval size (bp)");
        d.setYLab("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalA);
        d = new DensityPlot(BList.toArray());
        d.setSmoothN(5000);
        d.setTitle("Interval distribution of fd test in B subgenome");
        d.setXLab("Interval size (bp)");
        d.setYLab("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalB);

        d = new DensityPlot(DList.toArray());
        d.setSmoothN(5000);
        d.setTitle("Interval distribution of fd test in D subgenome");
        d.setXLab("Interval size (bp)");
        d.setYLab("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalD);
    }

}
