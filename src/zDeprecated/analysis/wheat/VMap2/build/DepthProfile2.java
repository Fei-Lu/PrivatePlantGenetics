package zDeprecated.analysis.wheat.VMap2.build;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import pgl.AppUtils;
import pgl.graph.r.DensityPlot;
import pgl.graph.r.ScatterPlot;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

class DepthProfile2 {

    public DepthProfile2 () {
//        this.mkTaxaBamDGenome();
//        this.mkPopDepInput();
//        this.mkPerlBatch();
//        this.modifyPopDepAB();
//        this.samplePopDep();
//        this.plotPopDep();



//        this.densityFilter();
        this.mkReliableSites();
//        this.profileReliableSites();
    }

    public void profileReliableSites () {
        String inDirS = "/Volumes/VMap2_Fei/reliableSites/round_02/ABD_intersect_80";
        String outfileS = "/Volumes/VMap2_Fei/reliableSites/round_02/reliable_count_80.txt";
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        int[] chromLengths = new int[fList.size()];
        int[] reliableNumber = new int[fList.size()];
        fList.parallelStream().forEach(f -> {
            int chr = Integer.parseInt(f.getName().split("_")[0].replace("chr", ""));
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                int length = 0;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    length++;
                    if (Integer.parseInt(temp) == 1) cnt++;
                }
                br.close();
                chromLengths[chr-1] = length;
                reliableNumber[chr-1] = cnt;
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName());
        });
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tChromLength\tReliableCount\tRatio");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < chromLengths.length; i++) {
                sb.setLength(0);
                sb.append(i+1).append("\t").append(chromLengths[i]).append("\t").append(reliableNumber[i]).
                        append("\t").append((float)((double)reliableNumber[i]/chromLengths[i]));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkReliableSites () {
        String abInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/AB/";
        String abd_abInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/ABD/";
        String abd_dInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/ABD/";
        String dInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/D/";

        String abInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/ab_popdep_sample.txt.gz";
        String abd_abInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_ab_popdep_sample.txt.gz";
        String abd_dInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_d_popdep_sample.txt.gz";
        String dInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/d_popdep_sample.txt.gz";

        String outDirS = "/Volumes/VMap2_Fei/reliableSites/round_02/ABD_intersect_80";
        new File(outDirS).mkdir();
        int binNum = 100;
        double proportion = 0.8;

        Grid abGrid = getGrid(2, 8, 2, 8, binNum, abInfileS);
        Grid abd_abGrid = getGrid(5, 14, 3, 8, binNum, abd_abInfileS);
        Grid abd_dGrid = getGrid(5, 14, 3, 8, binNum, abd_dInfileS);
        Grid dGrid = getGrid(7, 13, 3, 9, binNum, dInfileS);
        this.mkReliableSites(abGrid, abd_abGrid, abInDirS, abd_abInDirS, outDirS, proportion);
        this.mkReliableSites(dGrid, abd_dGrid, dInDirS, abd_dInDirS, outDirS, proportion);
    }

    private void mkReliableSites (Grid gr1, Grid gr2, String inDirS1, String inDirS2, String outDirS, double proportion) {
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS1, ".gz");
        List<File> f2List = IOUtils.getFileListInDirEndsWith(inDirS2, ".gz");
        fList.parallelStream().forEach(f -> {
            String header = f.getName().split("_")[0];
            String outfileS = new File(outDirS, header+"_"+String.valueOf(proportion)+"_reliable.txt.gz").getAbsolutePath();
            File f2 = null;
            for (int i = 0; i < f2List.size(); i++) {
                if (f2List.get(i).getName().startsWith(header)) {
                    f2 = f2List.get(i);
                    break;
                }
            }
            int indexThresh1 = gr1.getOrderIndexOfProportionOfSite(proportion);
            int indexThresh2 = gr2.getOrderIndexOfProportionOfSite(proportion);
            int v1 = 0;
            int v2 = 0;
            double x;
            double y;
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedReader br2 = IOUtils.getTextGzipReader(f2.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("IfReliable(0/1)");
                bw.newLine();
                String temp1 = br1.readLine();
                String temp2 = br2.readLine();
                List<String> l = new ArrayList<>();
                while ((temp1 = br1.readLine()) != null) {
                    temp2 = br2.readLine();
                    l = PStringUtils.fastSplit(temp1);
                    x = Double.parseDouble(l.get(1));
                    y = Double.parseDouble(l.get(2));
                    if (gr1.isHighDensity(x, y, indexThresh1)) v1 = 1;
                    else v1 = 0;
                    l = PStringUtils.fastSplit(temp2);
                    x = Double.parseDouble(l.get(1));
                    y = Double.parseDouble(l.get(2));
                    if (gr2.isHighDensity(x, y, indexThresh2)) v2 = 1;
                    else v2 = 0;
                    if (v1*v2 == 1) {
                        bw.write("1");
                    }
                    else {
                        bw.write("0");
                    }
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(outfileS);
        });
    }

    private Grid getGrid (double xlow, double xhigh, double ylow, double yhigh, int binNum, String infileS) {
        Grid gr = new Grid(xlow, xhigh, ylow, yhigh, binNum);
        RowTable<String> t = new RowTable<>(infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,2), t.getCellAsDouble(i,3));
        }
        gr.buildHashMap();
        return gr;
    }

    public void densityFilter () {
        double proportionOfSite = 0.80;
        int binNum = 100;
        String abInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/ab_popdep_sample.txt.gz";
        String abd_abInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_ab_popdep_sample.txt.gz";
        String abd_dInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_d_popdep_sample.txt.gz";
        String dInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/d_popdep_sample.txt.gz";
        String abOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/ab/scatter_p80.pdf";
        String abd_abOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/abd_ab/scatter_p80.pdf";
        String abd_dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/abd_d/scatter_p80.pdf";
        String dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/d/scatter_p80.pdf";
        this.densityFilter(abInfileS, abOutfileS, proportionOfSite, 2, 8, 2, 8, binNum, "AB");
        this.densityFilter(abd_abInfileS, abd_abOutfileS, proportionOfSite, 5, 14, 3, 8, binNum, "ABD_AB");
        this.densityFilter(abd_dInfileS, abd_dOutfileS, proportionOfSite, 5, 14, 3, 8, binNum, "ABD_D");
        this.densityFilter(dInfileS, dOutfileS, proportionOfSite, 7, 13, 3, 9, binNum, "D");
    }
    
    
    private void densityFilter (String infileS, String outfileS, double proportion, double xlow, double xhigh, double ylow, double yhigh, int binNum, String genomeType) {
        RowTable<String> t = new RowTable<>(infileS);
        Grid gr = new Grid(xlow, xhigh, ylow, yhigh, binNum);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,2), t.getCellAsDouble(i,3));
        }
        gr.buildHashMap();
        DoubleArrayList xList = new DoubleArrayList();
        DoubleArrayList yList = new DoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportion);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double x = t.getCellAsDouble(i, 2);
            double y = t.getCellAsDouble(i, 3);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
        }

        ScatterPlot s = new ScatterPlot(xList.toDoubleArray(), yList.toDoubleArray());
        s.setTitle(genomeType);
        s.setColor(255,0, 0, 1);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(outfileS);
    }

    public void plotPopDep () {
        String abInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/ab_popdep_sample.txt.gz";
        String abd_abInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_ab_popdep_sample.txt.gz";
        String abd_dInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_d_popdep_sample.txt.gz";
        String dInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/d_popdep_sample.txt.gz";
        String abOutDirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/ab";
        String abd_abOutDirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/abd_ab";
        String abd_dOutDirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/abd_d";
        String dOutDirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepPlot/dd";
        this.plotPopDep(abInfileS, abOutDirS, "AB");
        this.plotPopDep(abd_abInfileS, abd_abOutDirS, "ABD_AB");
        this.plotPopDep(abd_dInfileS, abd_dOutDirS, "ABD_D");
        this.plotPopDep(dInfileS, dOutDirS, "D");
    }

    private void plotPopDep (String infileS, String outDirS, String genomeType) {
        new File(outDirS).mkdir();
        DoubleArrayList meanList = new DoubleArrayList();
        DoubleArrayList sdList = new DoubleArrayList();
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                meanList.add(Double.parseDouble(l.get(2)));
                sdList.add(Double.parseDouble(l.get(3)));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        double[][] values = new double[2][];
        values[0] = meanList.toDoubleArray();
        values[1] = sdList.toDoubleArray();
        String meanOutfileS = new File (outDirS, "mean.pdf").getAbsolutePath();
        String sdOutFileS = new File (outDirS, "sd.pdf").getAbsolutePath();
        String scatterFileS = new File (outDirS, "scatter.pdf").getAbsolutePath();
        DensityPlot d = new DensityPlot(values[0]);
        d.setTitle(genomeType);
        d.setXLim(0, 20);
        d.setXLab("Mean of depth");
        d.setYLab("Density");
        d.setSmoothN(50000);
        d.saveGraph(meanOutfileS);

        d = new DensityPlot(values[1]);
        d.setTitle(genomeType);
        d.setXLim(0, 20);
        d.setXLab("SD of depth");
        d.setYLab("Density");
        d.setSmoothN(50000);
        d.saveGraph(sdOutFileS);

        ScatterPlot s = new ScatterPlot(values[0], values[1]);
        s.setTitle(genomeType);

        s.setXLim(0, 15);
        s.setYLim(0, 15);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.setColor(255, 0, 0, 1);
        s.saveGraph(scatterFileS);
    }

    public void samplePopDep () {
        int sampleSize = 100000;
        String abSource = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/AB";
        String abdSource = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/ABD";
        String dSource = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/D";
        String abOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/ab_popdep_sample.txt.gz";
        String abd_abOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_ab_popdep_sample.txt.gz";
        String abd_dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/abd_d_popdep_sample.txt.gz";
        String dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popDepSample/d_popdep_sample.txt.gz";
        int[] aChrIDs = RefV1Utils.getChrIDsOfSubgenomeA();
        int[] bChrIDs = RefV1Utils.getChrIDsOfSubgenomeB();
        int[] dChrIDs = RefV1Utils.getChrIDsOfSubgenomeD();
        int[] abChrIDs = new int[aChrIDs.length+bChrIDs.length];
        System.arraycopy(aChrIDs, 0, abChrIDs, 0, aChrIDs.length);
        System.arraycopy(bChrIDs, 0, abChrIDs, aChrIDs.length, bChrIDs.length);
        Arrays.sort(abChrIDs);
        this.samplePopDep(abSource, abOutfileS, abChrIDs, sampleSize);
        this.samplePopDep(abdSource, abd_abOutfileS, abChrIDs, sampleSize);
        this.samplePopDep(abdSource, abd_dOutfileS, dChrIDs, sampleSize);
        this.samplePopDep(dSource, dOutfileS, dChrIDs, sampleSize);
    }

    private void samplePopDep (String source, String outfileS, int[] chrIDs, int sampleSize) {
        Arrays.sort(chrIDs);
        List<File> fList = IOUtils.getFileListInDirEndsWith(source, ".gz");
        int[] siteNumbers = new int[chrIDs.length];
        int[] chrLengths = new int[chrIDs.length];
        double totalLength = 0;
        for (int i = 0; i < chrIDs.length; i++) {
            chrLengths[i] = RefV1Utils.getChrIDLength(chrIDs[i]);
            totalLength+=chrLengths[i];
        }
        for (int i = 0; i < chrIDs.length; i++) {
            siteNumbers[i] = (int)((chrLengths[i]/totalLength)*sampleSize);
        }
        List<String>[] resultLists = new ArrayList[chrIDs.length];
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < chrIDs.length; i++) {
            indexList.add(i);
            resultLists[i] = new ArrayList<String>();
        }
        indexList.parallelStream().forEach(i -> {
            String query = PStringUtils.getNDigitNumber(3, chrIDs[i]);
            String infileS = null;
            for (int j = 0; j < fList.size(); j++) {
                if (fList.get(j).getName().contains(query)) {
                    infileS = fList.get(j).getAbsolutePath();
                }
            }
            List<String> resultList = new ArrayList<>();
            int step = chrLengths[i]/siteNumbers[i];
            try {
                BufferedReader br = IOUtils.getTextGzipReader(infileS);
                String temp = br.readLine();
                int cnt = -1;
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%step != 0) continue;
                    sb.setLength(0);
                    sb.append(chrIDs[i]).append("\t").append(temp);
                    resultList.add(sb.toString());
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            resultLists[i] = resultList;
            System.out.println(chrIDs[i]);
        });
        try {
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            bw.write("Chr\tPos\tDepth_Mean\tDepth_SD");
            bw.newLine();
            for (int i = 0; i < resultLists.length; i++) {
                for (int j = 0; j < resultLists[i].size(); j++) {
                    bw.write(resultLists[i].get(j));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            System.out.println(outfileS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void modifyPopDepAB () {
        String sourceAB = "/Volumes/VMap2_Fei/popdep_vmap2/round_01/AB";
        String sourceABD = "/Volumes/VMap2_Fei/popdep_vmap2/round_01/ABD";
        String outAB = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/AB";
        String outABD = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/ABD";
        this.modifyPopDep(sourceAB, outAB);
        this.modifyPopDep(sourceABD, outABD);
    }

    private void modifyPopDep (String sourceDirS, String outDirS) {
        new File(outDirS).mkdir();
        List<File> fList = IOUtils.getFileListInDirEndsWith(sourceDirS, ".gz");
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                List<String> l = null;
                String temp = null;
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(temp);
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(outfileS);
        });
    }

    public void mkPerlBatch () {
//        java -Xmx100g -jar PlantGenetics.jar -a PopDep -p ./popdep_input/036_parameter_popdep.txt > log.txt
        String perlfile = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/run_popdep.pl";
        int[] chrIDs = RefV1Utils.getChrIDsOfSubgenomeD();
        StringBuilder sb = new StringBuilder();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(perlfile);
            for (int i = 0; i < chrIDs.length; i++) {
                sb.setLength(0);
                String currentChr = PStringUtils.getNDigitNumber(3, chrIDs[i]);
                String filename = currentChr + "_parameter_popdep.txt";
                String logname = currentChr + "_log.txt";
                sb.append("system(\"").append("java -Xmx100g -jar PlantGenetics.jar -a PopDep -p ./popdep_input/");
                sb.append(filename).append(" > ").append(logname).append("\");");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void mkPopDepInput () {
        String sampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/source/parameters_popdep.txt";
        String taxaBamFileS = "/data1/home/feilu/taxaBam.txt";
        String outDirS = "/data1/home/feilu/out";
        String samtoolsPath = "/data1/programs/samtools-1.8/samtools";
        String parametersDirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popdep_input";
        ArrayList<String> plinelist = new ArrayList<>();
        int[] chrIDs = RefV1Utils.getChrIDsOfSubgenomeD();
        for (int i = 0; i < chrIDs.length; i++) {
            plinelist.clear();
            plinelist.add(taxaBamFileS);
            plinelist.add(String.valueOf(chrIDs[i]));
            plinelist.add(String.valueOf(RefV1Utils.getChrIDLength(chrIDs[i])));
            plinelist.add(samtoolsPath);
            plinelist.add("32");
            String currentChr = PStringUtils.getNDigitNumber(3, chrIDs[i]);
            String resultFileS = new File (outDirS, currentChr+"_popdep.txt.gz").getAbsolutePath();
            String parameterFileS = new File (parametersDirS, currentChr+"_parameter_popdep.txt").getAbsolutePath();
            plinelist.add(resultFileS);
            AppUtils.creatParameterFile(sampleFileS, null, plinelist, parameterFileS);
        }

    }

    private void mkTaxaBamDGenome () {
        String dFileSYafei = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/source/yafei.txt";
        String dFileSXuebo = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/source/xuebo.txt";
        String dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popdep_input/taxaBam.txt";
        RowTable<String> t = new RowTable<>(dFileSYafei);
        List<String> l = t.getColumn(1);
        HashSet<String> s = new HashSet<>(l);
        t = new RowTable<>(dFileSXuebo);
        l = t.getColumn(5);
        s.addAll(l);
        l = new ArrayList<>(s);
        Collections.sort(l);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(dOutfileS);
            bw.write("Taxa	Bams(A list of bams of the taxon, seperated by the delimiter of Tab)");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < l.size(); i++) {
                sb.setLength(0);
                sb.append("D").append(i+1).append("\t").append(l.get(i));
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
}
