package analysis.wheat.VMap1;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.TIntHashSet;
import htsjdk.samtools.util.IOUtil;
import pgl.graph.r.DensityPlot;
import pgl.graph.r.Histogram;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
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
import java.util.*;

class Introgression {

    public Introgression () {
//        this.intervalSize();
//        this.findIndividualWithMaxFd();
//        this.maxFdHist();
        this.findIndividualWithMinIBSDistance();
    }

    public void findIndividualWithMinIBSDistance () {
        String vmap1DirS = "/Volumes/Fei_HDD_Mac/VMap1.0/VMapI";
        String maxFdDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/maxFd";
        String minIBSDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/minIBSD";
        String convertFileS = "/Users/feilu/Documents/analysisH/vmap1/fd/landrace_ID_convert.txt";
        RowTable<String> t = new RowTable<>(convertFileS);
        HashMap<Integer, String> landMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            landMap.put(i+1, t.getCell(i,2));
        }
        List<File> maxfdDirs = IOUtils.getDirListInDir(maxFdDirS);
        for (int i = 0; i < maxfdDirs.size(); i++) {
            String outDirS = new File (minIBSDirS, maxfdDirs.get(i).getName().split("_")[0]+"_minIBSD").getAbsolutePath();
            this.outputMinIBSDistance(maxfdDirs.get(i).getAbsolutePath(), outDirS, vmap1DirS, landMap);
        }
    }

    private void outputMinIBSDistance (String inDirS, String outDirS, String vcfDirS, HashMap<Integer, String> landMap) {
        File outDir = new File (outDirS);
        outDir.mkdir();
        List<File> infiles = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        Ranges rs = null;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infiles.get(0).getAbsolutePath());
            String temp = br.readLine();
            List<Range> rList = new ArrayList<>();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                Range r = new Range (Integer.parseInt(l.get(0)), Integer.parseInt(l.get(1)), Integer.parseInt(l.get(2)));
                rList.add(r);
            }
            rs = new Ranges(rList);
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        List<String> outfiles = new ArrayList<>();
        for (int i = 0; i < infiles.size(); i++) {
            outfiles.add(new File (outDirS, infiles.get(i).getName().split("_")[0]+"_minIBSD.txt").getAbsolutePath());
        }
        try {
            BufferedReader[] brs = new BufferedReader[infiles.size()];
            BufferedWriter[] bws = new BufferedWriter[outfiles.size()];
            int[] landID = new int[infiles.size()];
            for (int i = 0; i < brs.length; i++) {
                landID[i] = Integer.parseInt(infiles.get(i).getName().split("_")[0]);
                brs[i] = IOUtils.getTextGzipReader(infiles.get(i).getAbsolutePath());
                brs[i].readLine();
                bws[i] = IOUtils.getTextWriter(outfiles.get(i));
                bws[i].write("Chr\tStart\tEnd\tMinIBSDistance\tTaxa");
                bws[i].newLine();
            }
            String temp = null;
            List<String> l = new ArrayList<>();
            GenotypeGrid gt = null;
            StringBuilder sb = new StringBuilder();
            while ((temp = brs[0].readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int currentChr = Integer.parseInt(l.get(0));
                int currentStart = Integer.parseInt(l.get(1));
                int currentEnd = Integer.parseInt(l.get(2));
                if (gt == null || gt.getChromosome(0) != currentChr) {
                    String vcfFileS = new File (vcfDirS, "chr"+String.valueOf(currentChr)+".vcf.gz").getAbsolutePath();
                    gt = new GenotypeGrid(vcfFileS, GenoIOFormat.VCF_GZ);
                    gt.sortByTaxa();
                }
                int startIndex = gt.getSiteIndex((short)currentChr, currentStart);
//                if (startIndex < 0) startIndex = -startIndex-1;
                int endIndex = gt.getSiteIndex((short)currentChr, currentEnd);
//                if (endIndex < 0) endIndex = -endIndex-1;
                for (int i = 0; i < brs.length; i++) {
                    if (i != 0) temp = brs[i].readLine();
                    this.processIBSOutputLine(currentChr, currentStart, currentEnd, temp, gt, landMap.get(landID[i]), sb);
                    bws[i].write(sb.toString());
                    bws[i].newLine();
                }
            }
            for (int i = 0; i < brs.length; i++) {
                bws[i].flush();
                bws[i].close();
                brs[i].close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void processIBSOutputLine (int chr, int start, int end, String inputLine, GenotypeGrid gt, String currentTaxon, StringBuilder sb) {
        sb.setLength(0);
        sb.append(chr).append("\t").append(start).append("\t").append(end).append("\t");
        List<String> l = PStringUtils.fastSplit(inputLine);
        float fd = Float.parseFloat(l.get(3));
        if (fd == 1.0) {
            int n = l.size()-4;
            float[] ds = new float[n];
            Arrays.fill(ds, 1);
            int currentIndex = gt.getTaxonIndex(currentTaxon);
            for (int i = 0; i < n; i++) {
                int nextIndex = gt.getTaxonIndex(l.get(i+4));
                ds[i] = gt.getIBSDistance(currentIndex, nextIndex, start, end);
            }
            float minIBSD = 1;
            for (int i = 0; i < n; i++) {
                if (ds[i] < minIBSD) minIBSD = ds[i];
            }
            sb.append(minIBSD);
            for (int i = 0; i < n; i++) {
                if (ds[i] == minIBSD) {
                    sb.append("\t").append(l.get(i+4));
                }
            }
        }
        else {
            sb.append("NA");
        }
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
