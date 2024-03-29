package zDeprecated.analysis.wheat.VMap1;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.graph.r.BoxPlot;
import pgl.graph.r.DensityPlot;
import pgl.graph.r.Histogram;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeTable;
import pgl.infra.range.Range;
import pgl.infra.range.Ranges;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Dyad;
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
//        this.findIndividualWithMinIBSDistance();
//        this.calculateProportionOfIntrogressionOfIndividual();
//        this.calculateContributionToIndividualIntrogression();

//        this.plotIBSDistanceDistribution();
//        this.plotIBSDistanceByRegion();

    }

    public void plotIBSDistanceByRegion () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/minIBSD_dis";
        String[] regions = {"EU", "WA", "EA"};
        String[] subgenomes = {"A", "B", "D"};
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < subgenomes.length; i++) {
            double[][] values = new double[regions.length][];
            for (int j = 0; j < regions.length; j++) {
                sb.setLength(0);
                sb.append(regions[j]).append("_").append(subgenomes[i]).append("_IBSD_dis.txt");
                String infileS = new File (inDirS, sb.toString()).getAbsolutePath();
                RowTable<String> t = new RowTable<>(infileS);
                values[j] = t.getColumnAsDoubleArray(0);
                values[j] = this.logTransformIBSD(values[j]);
            }
            sb.setLength(0);
            sb.append(subgenomes[i]).append("_IBSD_box.pdf");
            String outfileS = new File (inDirS, sb.toString()).getAbsolutePath();
            BoxPlot b = new BoxPlot(values, regions);
            b.setTitle("IBSD_"+subgenomes[i]);
            b.setYLabel("IBS Distance");
            b.saveGraph(outfileS);
        }
    }

    private double[] logTransformIBSD (double[] xs) {
        double[] values = new double[xs.length];
        for (int i = 0; i < values.length; i++) {
            values[i] = Math.log10(xs[i]+0.00001);
        }
        return values;
    }

    public void plotIBSDistanceDistribution () {
        String landraceRegionFileS = "/Users/feilu/Documents/analysisH/vmap1/fd/Landrace_region.txt";
        String minIBSDDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/minIBSD";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/minIBSD_dis";
        float fdThresh = (float)0.5;
        RowTable<String> t = new RowTable<>(landraceRegionFileS);
        int[] euIDs = this.getRegionIDs(t, "EU");
        int[] waIDs = this.getRegionIDs(t, "WA");
        int[] eaIDs = this.getRegionIDs(t, "EA");
        double[][] ibsA = this.getIBSDofTaxon(minIBSDDirS, "A", fdThresh);
        double[][] ibsB = this.getIBSDofTaxon(minIBSDDirS, "B", fdThresh);
        double[][] ibsD = this.getIBSDofTaxon(minIBSDDirS, "D", fdThresh);
        int[][] IDs = new int[3][];
        IDs[0] = euIDs;IDs[1] = waIDs;IDs[2] = eaIDs;
        double[][][] ibs = new double[3][][];
        ibs[0] = ibsA;ibs[1] = ibsB;ibs[2] = ibsD;
        String[] regions = {"EU", "WA", "EA"};
        String[] subgenomes = {"A", "B", "D"};
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < regions.length; i++) {
            for (int j = 0; j < subgenomes.length; j++) {
                sb.setLength(0);
                sb.append(regions[i]).append("_").append(subgenomes[j]).append("_IBSD_dis.pdf");
                String outfileS = new File (outDirS, sb.toString()).getAbsolutePath();
                String outtextFileS = outfileS.replaceFirst(".pdf", ".txt");
                TDoubleArrayList dList = new TDoubleArrayList();
                for (int k = 0; k < IDs[i].length; k++) {
                    int taxonIndex = IDs[i][k]-1;
                    dList.addAll(ibs[j][taxonIndex]);
                }
                double[] diss = dList.toArray();
                Histogram h = new Histogram(diss);
                h.setTitle("IBS distance of " + regions[i] + "_" + subgenomes[j]);
                h.setXLim(0, 0.3);
                h.setXLabel("IBS distance");
                h.setYLabel("Proportion");
                h.setBreakNumber(50);
                h.saveGraph(outfileS);
                try {
                    BufferedWriter bw = IOUtils.getTextWriter(outtextFileS);
                    bw.write("IBSDistance");
                    bw.newLine();
                    for (int k = 0; k < diss.length; k++) {
                        bw.write(String.valueOf((float)diss[k]));
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
    }

    private double[][] getIBSDofTaxon (String minIBSDDirS,  String subgenome, float fdThresh) {
        List<File> ibsDirs = IOUtils.getDirListInDir(minIBSDDirS);
        double[][] value = null;
        for (int i = 0; i < ibsDirs.size(); i++) {
            if (ibsDirs.get(i).getName().split("_")[0].equals(subgenome)) {
                List<File> fList = IOUtils.getFileListInDirEndsWith(ibsDirs.get(i).getAbsolutePath(), ".txt");
                value = new double[fList.size()][];
                for (int j = 0; j < fList.size(); j++) {
                    int taxonIndex = Integer.parseInt(fList.get(j).getName().split("_")[0])-1;
                    TDoubleArrayList dList = new TDoubleArrayList();
                    try {
                        BufferedReader br = IOUtils.getTextReader(fList.get(j).getAbsolutePath());
                        String temp = br.readLine();
                        List<String> l = new ArrayList<>();
                        while ((temp = br.readLine()) != null) {
                            l = PStringUtils.fastSplit(temp);
                            if (l.get(4).startsWith("N")) continue;
                            float fd = Float.parseFloat(l.get(3));
                            if (fd < fdThresh) continue;
                            dList.add(Double.parseDouble(l.get(4)));
                        }
                        br.close();
                        value[taxonIndex] = dList.toArray();
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }
        }
        return value;
    }

    private int[] getRegionIDs (RowTable<String> t, String region) {
        TIntArrayList idList = new TIntArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCell(i, 2).equals(region)) {
                idList.add(Integer.parseInt(t.getCell(i,0)));
            }
        }
        int[] ids = idList.toArray();
        Arrays.sort(ids);
        return ids;
    }

    public void calculateContributionToIndividualIntrogression() {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/minIBSD";
        String taxaFileS = "/Users/feilu/Documents/analysisH/vmap1/fd/taxa_info.txt";
        String landFileS = "/Users/feilu/Documents/analysisH/vmap1/fd/Landrace_region.txt";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap1/fd/individual_intro/individual_contribution.txt";
        RowTable<String> t = new RowTable<>(taxaFileS);
        String[] types = {"Urartu", "Wild_emmer", "Domesticated_emmer", "Strangulata", "Free_threshing"};
        Arrays.sort(types);
        HashMap<String, String> taxaTypeMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(types, t.getCell(i,9));
            if (index < 0) continue;
            taxaTypeMap.put(t.getCell(i, 0), types[index]);
        }
        t = new RowTable<>(landFileS);
        String[] subgenomes = {"A", "B", "D"};
        List<File> dirs = new ArrayList<>();
        for (int i = 0; i < subgenomes.length; i++) {
            dirs.add(new File (inDirS, subgenomes[i]+"_minIBSD"));
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = "TaxaID\tTaxa\tRegion\tUrautu\tWE\tDE\tFree\tStrangulata";
            StringBuilder sb = new StringBuilder("TaxaID\tTaxa\tRegion");
            for (int i = 0; i < types.length; i++) {
                sb.append("\t").append(types[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append(t.getCell(i, 0)).append("\t");
                sb.append(t.getCell(i, 1)).append("\t");
                sb.append(t.getCell(i, 2));
                double[] lengths = new double[types.length];
                long total = 0;
                String temp = null;
                for (int j = 0; j < dirs.size(); j++) {
                    String infileS = new File (dirs.get(j).getAbsolutePath(), String.valueOf(i+1)+"_minIBSD.txt").getAbsolutePath();
                    try {
                        BufferedReader  br = IOUtils.getTextReader(infileS);
                        temp = br.readLine();
                        List<String> l = new ArrayList<>();
                        while ((temp = br.readLine()) != null) {
                            l = PStringUtils.fastSplit(temp);
                            int len = Integer.parseInt(l.get(2))- Integer.parseInt(l.get(1));
                            total+=len;
                            if (l.get(4).startsWith("N")) continue;
                            double fd = Double.parseDouble(l.get(3));
                            int n = l.size()-5;
                            boolean[] ifOut = new boolean[types.length];
                            for (int k = 0; k < n; k++) {
                                String query = taxaTypeMap.get(l.get(k+5));
                                int index = Arrays.binarySearch(types, query);
                                if (!ifOut[index]) {
                                    ifOut[index] = true;
                                    lengths[index]+=fd*len;
                                }
                            }
                        }
                        br.close();
                    }
                    catch (Exception e) {
//                        System.out.println(infileS);
//                        System.out.println(temp);
                        e.printStackTrace();
                        //System.exit(1);
                    }
                }
                for (int j = 0; j < types.length; j++) {
                    sb.append("\t").append((float)(lengths[j]/total));
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

    public void calculateProportionOfIntrogressionOfIndividual() {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap1/fd/minIBSD";
        String landRegionFileS = "/Users/feilu/Documents/analysisH/vmap1/fd/Landrace_region.txt";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap1/fd/individual_intro/individual_proportion.txt";
        String[] subgenomes = {"A", "B", "D"};
        List<File> dirs = new ArrayList<>();
        RowTable<String> t = new RowTable<>(landRegionFileS);
        Dyad<Long, float[]>[] results = new Dyad[subgenomes.length];
        for (int i = 0; i < subgenomes.length; i++) {
            dirs.add(new File (inDirS, subgenomes[i]+"_minIBSD"));
            results[i] = this.getProportionOfSubgenome(dirs.get(i), t.getRowNumber());
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("TaxaID\tTaxa\tRegion\tA_Introgression\tB_introgression\tD_introgression\tAB_introgression\tAll_introgression");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append(t.getCell(i, 0)).append("\t");
                sb.append(t.getCell(i, 1)).append("\t");
                sb.append(t.getCell(i, 2)).append("\t");
                long total = 0;
                long abTotal = 0;
                for (int j = 0; j < subgenomes.length; j++) {
                    sb.append(results[j].getSecondElement()[i]).append("\t");
                    total+=results[j].getFirstElement();
                    if (j== 2) continue;
                    abTotal+=results[j].getFirstElement();
                }
                double pro = 0;
                double abPro = 0;
                for (int j = 0; j < subgenomes.length; j++) {
                    pro+=results[j].getSecondElement()[i]*((double)results[j].getFirstElement()/total);
                    if (j== 2) continue;
                    abPro+=results[j].getSecondElement()[i]*((double)results[j].getFirstElement()/abTotal);
                }
                sb.append((float)abPro).append("\t").append((float)pro);
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

    private Dyad<Long, float[]> getProportionOfSubgenome (File dir, int taxaNum) {
        long total = 0;
        TIntArrayList regionList = new TIntArrayList();
        File[] fs = new File[taxaNum];
        float[] values = new float[taxaNum];
        for (int i = 0; i < fs.length; i++) {
            fs[i] = new File (dir, String.valueOf(i+1)+"_minIBSD.txt");
            total = 0;
            double actual = 0;
            try {
                BufferedReader br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = br.readLine();
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int length = Integer.parseInt(l.get(2))- Integer.parseInt(l.get(1));
                    total+=length;
                    if (l.get(4).startsWith("N")) continue;
                    double fd = Double.parseDouble(l.get(3));
                    actual+=(fd*length);
                }
                values[i] = (float)((double)actual/total);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        return new Dyad<>(total, values);
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
        new File(minIBSDirS).mkdir();
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
                bws[i].write("Chr\tStart\tEnd\tMaxFd\tMinIBSDistance\tTaxa");
                bws[i].newLine();
            }
            String temp = null;
            List<String> l = new ArrayList<>();
            GenotypeTable gt = null;
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

    private void processIBSOutputLine (int chr, int start, int end, String inputLine, GenotypeTable gt, String currentTaxon, StringBuilder sb) {
        int startIndex = gt.getSiteIndex((short)chr, start);
//                if (startIndex < 0) startIndex = -startIndex-1;
        int endIndex = gt.getSiteIndex((short)chr, end);
//                if (endIndex < 0) endIndex = -endIndex-1;
        sb.setLength(0);
        sb.append(chr).append("\t").append(start).append("\t").append(end).append("\t");
        List<String> l = PStringUtils.fastSplit(inputLine);
        sb.append(l.get(3)).append("\t");
        float fd = Float.parseFloat(l.get(3));
        //set the threshold of fd here
        if (fd > 0.5) {
            int n = l.size()-4;
            float[] ds = new float[n];
            Arrays.fill(ds, 1);
            int currentIndex = gt.getTaxonIndex(currentTaxon);
            for (int i = 0; i < n; i++) {
                int nextIndex = gt.getTaxonIndex(l.get(i+4));
                ds[i] = gt.getIBSDistance(currentIndex, nextIndex, startIndex, endIndex);
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
            h.setXLabel("fd");
            h.setYLabel("Proportion");
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
        d.setXLabel("Interval size (bp)");
        d.setYLabel("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalA);
        d = new DensityPlot(BList.toArray());
        d.setSmoothN(5000);
        d.setTitle("Interval distribution of fd test in B subgenome");
        d.setXLabel("Interval size (bp)");
        d.setYLabel("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalB);

        d = new DensityPlot(DList.toArray());
        d.setSmoothN(5000);
        d.setTitle("Interval distribution of fd test in D subgenome");
        d.setXLabel("Interval size (bp)");
        d.setYLabel("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalD);
    }

}
