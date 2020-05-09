package analysis.wheat.VMap2.build;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.app.popdep.PopDep;
import pgl.graphcis.r.DensityPlot;
import pgl.graphcis.r.ScatterPlot;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

class DepthProfile {

    public DepthProfile () {
//        this.mkTaxaBamFile();
//        this.mkChrLengthFile();
//        this.runStep1();
//        this.mkStep2FileS();
//        this.runStep2();
        //this.popdepPlotSample();
        //this.popdepPlot();
        //this.vcfPlotSample();
        //this.vcfPlot();
        //this.findPopDepMode();
//        this.densityFilter();
        //this.mkReliableGenotypeSite();
        //this.mkReliableIntersection();
        //this.intersectionCheck();
    }

    public void intersectionCheck () {
        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_sample.txt";

        String abFile1S = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abPopDep_sample_chr001.txt";
        String abFile2S = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abPopDep_sample_chr002.txt";
        String abdFile1S = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abdPopDep_sample_chr001.txt";
        String abdFile2S = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abdPopDep_sample_chr002.txt";
        String dFile1S = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/dPopDep_sample_chr005.txt";
        String dFile2S = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/dPopDep_sample_chr006.txt";

//        this.split(abSampleFileS, abFile1S, abFile2S);
//        this.split(abdSampleFileS, abdFile1S, abdFile2S);
//        this.split(dSampleFileS, dFile1S, dFile2S);

        String reliableChr001 = "/Volumes/VMap2_Fei/reliableSites/ABD_intersect/chr001_intersect_reliable.txt.gz";
        String reliableChr002 = "/Volumes/VMap2_Fei/reliableSites/ABD_intersect/chr002_intersect_reliable.txt.gz";
        String reliableChr005 = "/Volumes/VMap2_Fei/reliableSites/ABD_intersect/chr005_intersect_reliable.txt.gz";
        String reliableChr006 = "/Volumes/VMap2_Fei/reliableSites/ABD_intersect/chr006_intersect_reliable.txt.gz";

        String abDepthDensity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abPopDep_intersect_depth_density.pdf";
        String abScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abPopDep_intersect_scatter.pdf";
        String abdDepthDensity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abdPopDep_intersect_depth_density.pdf";
        String abdScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/abdPopDep_intersect_scatter.pdf";
        String dDepthDensity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/dPopDep_intersect_depth_density.pdf";
        String dScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/intersection/dPopDep_intersect_scatter.pdf";
        this.plotIntersection(abFile1S, abFile2S, reliableChr001, reliableChr002, abDepthDensity, abScatter, "AB");
        this.plotIntersection(abdFile1S, abdFile2S, reliableChr001, reliableChr002, abdDepthDensity, abdScatter, "ABD");
        this.plotIntersection(dFile1S, dFile2S, reliableChr005, reliableChr006, dDepthDensity, dScatter, "D");
    }

    private void plotIntersection (String file1, String file2, String reliableFile1, String reliableFile2, String densityPlot, String scatterPlot, String genomeType) {
        TDoubleArrayList depthList = new TDoubleArrayList();
        TDoubleArrayList sdList = new TDoubleArrayList();
        int totalCnt = 0;
        int depthFilterCnt = 0;
        int[] values =this.plotIntersection2(file1, reliableFile1, depthList, sdList);
        totalCnt+=values[0]; depthFilterCnt+=values[1];
        values = this.plotIntersection2(file2, reliableFile2, depthList, sdList);
        totalCnt+=values[0]; depthFilterCnt+=values[1];
        double[] x = depthList.toArray();
        double[] y = sdList.toArray();
        ScatterPlot s = new ScatterPlot(x, y);
        DensityPlot d = new DensityPlot(x);
        if (genomeType.equals("AB")) {
            s.setTitle("AB_intersect");
            s.setXLim(0, 8);
            s.setYLim(0, 8);
            s.setColor(255, 0, 0, 5);
            s.setXLab("Mean of depth");
            s.setYLab("SD of depth");
            s.saveGraph(scatterPlot);
            d.setXLim(0, 8);
            d.setXLab("Mean of depth");
            d.setYLab("Density");
            d.setTitle("AB_intersect");
            d.saveGraph(densityPlot);
        }
        else if (genomeType.equals("ABD")) {
            s.setTitle("ABD_intersect");
            s.setXLim(0, 20);
            s.setYLim(0, 12);
            s.setColor(255, 0, 0, 5);
            s.setXLab("Mean of depth");
            s.setYLab("SD of depth");
            s.saveGraph(scatterPlot);
            d.setXLim(0, 20);
            d.setXLab("Mean of depth");
            d.setYLab("Density");
            d.setTitle("ABD_intersect");
            d.saveGraph(densityPlot);
        }
        else if (genomeType.equals("D")) {
            s.setTitle("D_intersect");
            s.setXLim(0, 20);
            s.setYLim(0, 12);
            s.setColor(255, 0, 0, 5);
            s.setXLab("Mean of depth");
            s.setYLab("SD of depth");
            s.saveGraph(scatterPlot);
            d.setXLim(0, 20);
            d.setXLab("Mean of depth");
            d.setYLab("Density");
            d.setTitle("D_intersect");
            d.saveGraph(densityPlot);
        }
        System.out.println(genomeType);
        System.out.println((double)x.length/totalCnt);

    }

    private int[] plotIntersection2 (String file, String reliableFile, TDoubleArrayList depthList, TDoubleArrayList sdList) {
        RowTable<String> t = new RowTable<>(file);
        int[] positions = t.getColumnAsIntArray(0);
        boolean[] ifSelected = new boolean[positions.length];
        int depthFilterCnt = 0;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(reliableFile);
            String temp = br.readLine();
            int currentPos = 0;
            int value;
            int index;

            while ((temp = br.readLine()) != null) {
                currentPos++;
                if (currentPos%10000000 == 0) System.out.println(currentPos);
                value = Integer.parseInt(temp);
                if (value == 0) {
                    continue;
                }
                index = Arrays.binarySearch(positions, currentPos);
                if (index < 0) {

                    continue;
                }
                ifSelected[index] = true;
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!ifSelected[i]) continue;
            depthList.add(Double.parseDouble(t.getCell(i,1)));
            sdList.add(Double.parseDouble(t.getCell(i,2)));
        }
        int[] values = new int[2];
        values[0] = t.getRowNumber();
        values[1] = t.getRowNumber()-depthFilterCnt;
        return values;
    }

    private void split (String sourceFileS, String file1, String file2) {
        try {
            BufferedReader br = IOUtils.getTextReader(sourceFileS);
            String header  = br.readLine();
            String temp = null;
            BufferedWriter bw1 = IOUtils.getTextWriter(file1);
            BufferedWriter bw2 = IOUtils.getTextWriter(file2);
            bw1.write(header);
            bw1.newLine();
            bw2.write(header);
            bw2.newLine();
            List<String> l = new ArrayList<>();
            int current = -1;
            int next = -1;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                next = Integer.parseInt(l.get(0));
                if (next < current) {
                    bw2.write(temp);
                    bw2.newLine();
                    while ((temp = br.readLine()) != null) {
                        bw2.write(temp);
                        bw2.newLine();
                    }
                }
                else {
                    bw1.write(temp);
                    bw1.newLine();
                }
                current = next;
            }
            bw1.flush();
            bw1.close();
            bw2.flush();
            bw2.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkReliableIntersection () {
        String abDirS = "/Volumes/VMap2_Fei/reliableSites/AB/";
        String abdDirS = "/Volumes/VMap2_Fei/reliableSites/ABD/";
        String dDirS = "/Volumes/VMap2_Fei/reliableSites/D/";

        String resultDirS = "/Volumes/VMap2_Fei/reliableSites/ABD_intersect";

        this.mkReliable(abDirS, abdDirS, resultDirS);
        this.mkReliable(dDirS, abdDirS, resultDirS);

    }

    private void mkReliable (String source1DirS, String source2DirS, String resultDirS) {
        List<File> fList = IOUtils.getFileListInDirEndsWith(source1DirS, ".gz");
        Collections.sort(fList);
        fList.parallelStream().forEach(f -> {
            String source2FileS = new File (source2DirS, f.getName().split("_")[0]+"_ABD_reliable.txt.gz").getAbsolutePath();
            String outputFileS = new File (resultDirS, f.getName().split("_")[0]+"_intersect_reliable.txt.gz").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedReader br2 = IOUtils.getTextGzipReader(source2FileS);
                BufferedWriter bw = IOUtils.getTextGzipWriter(outputFileS);
                String temp1 = br1.readLine();
                String temp2 = br2.readLine();
                bw.write(temp1);
                bw.newLine();
                int cnt = 0;
                while ((temp1 = br1.readLine()) != null) {
                    temp2 = br2.readLine();
                    bw.write(String.valueOf(Integer.parseInt(temp1)*Integer.parseInt(temp2)));
                    bw.newLine();
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(cnt);
                }
                bw.flush();
                bw.close();
                br1.close();
                br2.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void mkReliableGenotypeSite () {
        String abInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/AB/";
        String abdInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/ABD/";
        String dInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/D/";

        String abOutDirS = "/Volumes/VMap2_Fei/reliableSites/AB";
        String abdOutDirS = "/Volumes/VMap2_Fei/reliableSites/ABD";
        String dOutDirS = "/Volumes/VMap2_Fei/reliableSites/D";

        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_sample.txt";

        double depthStart = 2;
        double depthEnd = 8;
        double SDStart = 2;
        double SDEnd = 8;
        int nBin = 100;
        double proportionOfSite = 0.7;

        this.mkReliable(abInDirS, abOutDirS, abSampleFileS, depthStart, depthEnd, SDStart, SDEnd, nBin, proportionOfSite);

        depthStart = 3;
        depthEnd = 15;
        SDStart = 3;
        SDEnd = 8;
        this.mkReliable(abdInDirS, abdOutDirS, abdSampleFileS, depthStart, depthEnd, SDStart, SDEnd, nBin, proportionOfSite);

        depthStart = 3;
        depthEnd = 17;
        SDStart = 3;
        SDEnd = 10;
        this.mkReliable(dInDirS, dOutDirS, dSampleFileS, depthStart, depthEnd, SDStart, SDEnd, nBin, proportionOfSite);

    }

    private void mkReliable (String inDirS, String outDirS, String sampleFileS, double depthStart, double depthEnd, double SDStart, double SDEnd, int nBin, double proportionOfSite) {
        RowTable<String> t = new RowTable<>(sampleFileS);
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS,".gz");
        fList.stream().forEach(f -> {
            String outfileS = f.getName().replaceFirst("popdep_vmap2.txt.gz", "reliable.txt.gz");
            outfileS = new File(outDirS, outfileS).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("IfReliable(0/1)");
                bw.newLine();
                List<String> l = new ArrayList<>();
                String temp = br.readLine();
                double x;
                double y;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    x = Double.parseDouble(l.get(1));
                    y = Double.parseDouble(l.get(2));
                    if (gr.isHighDensity(x, y, indexThresh)) {
                        bw.write("1");
                    }
                    else {
                        bw.write("0");
                    }
                    bw.newLine();
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(cnt);
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void densityFilter () {
        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_sample.txt";
        double proportionOfSite = 0.70;
        this.densityFilterPlotAB(proportionOfSite, abSampleFileS);
        this.densityFilterPlotABD(proportionOfSite, abdSampleFileS);
        this.densityFilterPlotD(proportionOfSite, dSampleFileS);

    }

    public void densityFilterPlotD(double proportionOfSite, String inputFileS) {
        String beforeFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/d_before_scatter.pdf";
        String afterFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/d_after_scatter.pdf";
        RowTable<String> t = new RowTable<>(inputFileS);
        double depthStart = 3;
        double depthEnd = 17;
        double SDStart = 3;
        double SDEnd = 10;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double x = t.getCellAsDouble(i, 1);
            double y = t.getCellAsDouble(i, 2);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
        }

        ScatterPlot s = new ScatterPlot(t.getColumnAsDoubleArray(1), t.getColumnAsDoubleArray(2));
        s.setTitle("D_before");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(beforeFileS);
        s = new ScatterPlot(xList.toArray(), yList.toArray());
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(afterFileS);
    }
    public void densityFilterPlotABD(double proportionOfSite, String inputFileS) {
        String beforeFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/abd_before_scatter.pdf";
        String afterFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/abd_after_scatter.pdf";
        RowTable<String> t = new RowTable<>(inputFileS);
        double depthStart = 3;
        double depthEnd = 15;
        double SDStart = 3;
        double SDEnd = 8;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double x = t.getCellAsDouble(i, 1);
            double y = t.getCellAsDouble(i, 2);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
        }

        ScatterPlot s = new ScatterPlot(t.getColumnAsDoubleArray(1), t.getColumnAsDoubleArray(2));
        s.setTitle("ABD_before");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(beforeFileS);
        s = new ScatterPlot(xList.toArray(), yList.toArray());
        s.setTitle("ABD_after");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.saveGraph(afterFileS);
    }

    public void densityFilterPlotAB(double proportionOfSite, String inputFileS) {
        String beforeFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/ab_before_scatter.pdf";
        String afterFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/ab_after_scatter.pdf";
        RowTable<String> t = new RowTable<>(inputFileS);
        double depthStart = 2;
        double depthEnd = 8;
        double SDStart = 2;
        double SDEnd = 8;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double x = t.getCellAsDouble(i, 1);
            double y = t.getCellAsDouble(i, 2);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
        }
        ScatterPlot s = new ScatterPlot(t.getColumnAsDoubleArray(1), t.getColumnAsDoubleArray(2));
        s.setTitle("AB_before");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 8);
        s.setYLim(0, 8);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(beforeFileS);
        s = new ScatterPlot(xList.toArray(), yList.toArray());
        s.setTitle("AB_after");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 8);
        s.setYLim(0, 8);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(afterFileS);
    }

    /**
     * By aoyue
     */
    public void findPopDepMode () {
        double abDepthMode = 5.144204;
        double abdDepthMode = 9.607552;
        double dDepthMode = 11.91677;
        double abSDMode = 4.226852;
        double abdSDMode = 5.162204;
        double dSDMode = 7.056306;
    }

    public void vcfPlot () {
        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abVCF_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdVCF_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dVCF_sample.txt";

        String abDepthDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abVCF_depth_density.pdf";
        String abDepthScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abVCF_depth_scatter.pdf";

        String abdDepthDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdVCF_depth_density.pdf";
        String abdDepthScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdVCF_depth_scatter.pdf";

        String dDepthDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dVCF_depth_density.pdf";
        String dDepthScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dVCF_depth_scatter.pdf";
        this.plotVCF(abSampleFileS, abDepthDesity, abDepthScatter, "AB");
        this.plotVCF(abdSampleFileS, abdDepthDesity, abdDepthScatter, "ABD");
        this.plotVCF(dSampleFileS, dDepthDesity, dDepthScatter, "D");

    }

    private void plotVCF (String infileS, String depthDensity, String depthScatter, String genomeType) {
        RowTable<String> t = new RowTable<>(infileS);
        double[] depth = t.getColumnAsDoubleArray(1);
        double[] depthSD = t.getColumnAsDoubleArray(2);

        DensityPlot d = new DensityPlot(depth);
        d.setTitle(genomeType);
        d.setXLim(0, 20);
        d.setXLab("Mean of depth");
        d.setYLab("Density");
        d.saveGraph(depthDensity);
        ScatterPlot s = new ScatterPlot(depth, depthSD);
        s.setTitle(genomeType);
        s.setXLim(0, 20);
        s.setYLim(0, 20);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.setColor(255, 0, 0, 20);
        s.saveGraph(depthScatter);
    }

    public void vcfPlotSample () {
        String abFile1 = "/Volumes/VMap2_Fei/vcf/001_fastcall/ab/chr001.ABgenome.vcf.gz";
        String abFile2 = "/Volumes/VMap2_Fei/vcf/001_fastcall/ab/chr002.ABgenome.vcf.gz";
        String abdFile1 = "/Volumes/VMap2_Fei/vcf/001_fastcall/abd/chr001.ABDgenome.vcf.gz";
        String abdFile2 = "/Volumes/VMap2_Fei/vcf/001_fastcall/abd/chr001.ABDgenome.vcf.gz";
        String dFile1 = "/Volumes/VMap2_Fei/vcf/001_fastcall/d/chr005.Dgenome.vcf.gz";
        String dFile2 = "/Volumes/VMap2_Fei/vcf/001_fastcall/d/chr006.Dgenome.vcf.gz";
        int sampleSize = 10000;
        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abVCF_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdVCF_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dVCF_sample.txt";
        this.sampleVCF(abFile1, abFile2, abSampleFileS, sampleSize);
        this.sampleVCF(abdFile1, abdFile2, abdSampleFileS, sampleSize);
        this.sampleVCF(dFile1, dFile2, dSampleFileS, sampleSize);
    }

    public void sampleVCF (String file1, String file2, String outfileS, int sampleSize) {
        File[] fs = new File[2];
        fs[0] = new File(file1);
        fs[1] = new File(file2);
        int[] chrs = new int[fs.length];
        int[] chrLengths = new int[fs.length];
        int[] siteCounts = new int[fs.length];
        int totalCount = 0;
        for (int i = 0; i < fs.length; i++) {
            chrs[i] = Integer.parseInt(fs[i].getName().split("\\.")[0].replaceFirst("chr", ""));
            chrLengths[i] = RefV1Utils.getChrIDLength(chrs[i]);
            try {
                BufferedReader br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) continue;
                    break;
                }
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(cnt);
                }
                siteCounts[i] = cnt;
                br.close();
                totalCount+=cnt;
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        int[] sampleCounts = new int[fs.length];
        for (int i = 0; i < fs.length; i++) {
            sampleCounts[i] = (int)Math.round((double)siteCounts[i]/totalCount*sampleSize);
            System.out.println(sampleCounts[i]);
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Position\tDepth_Mean\tDepth_SD");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                int step = siteCounts[i]/sampleCounts[i];
                BufferedReader br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) continue;
                    break;
                }
                int cnt = -1;
                List<String> l = new ArrayList<>();
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(cnt);
                    if (cnt%step == 0) {
                        TDoubleArrayList dList = new TDoubleArrayList();
                        l = PStringUtils.fastSplit(temp);
                        int n = l.size() - 9;
                        String current;
                        for (int j = 0; j < n; j++) {
                            current = l.get(j+9);
                            if (current.startsWith(".")) {
                                dList.add(0);
                            }
                            else {
                                String[] tem = current.split(":")[1].split(",");
                                int dep = 0;
                                for (int k = 0; k < tem.length; k++) {
                                    dep+=Integer.parseInt(tem[k]);
                                }
                                dList.add(dep);
                            }
                        }
                        sb.setLength(0);
                        DescriptiveStatistics d = new DescriptiveStatistics(dList.toArray());
                        sb.append(l.get(1)).append("\t").append((float)d.getMean()).append("\t").append((float)d.getStandardDeviation());
                        bw.write(sb.toString());
                        bw.newLine();

                    }
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("done");
    }

    public void popdepPlot () {
        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_sample.txt";

        String abDepthDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_depth_density.pdf";
        String abDepthSDDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_SD_density.pdf";
        String abDepthStanDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_depth_stan_density.pdf";
        String abDepthScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_depth_scatter.pdf";
        String abDepthStanScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_depth_stan_scatter.pdf";

        String abdDepthDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_depth_density.pdf";
        String abdDepthSDDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_SD_density.pdf";
        String abdDepthStanDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_depth_stan_density.pdf";
        String abdDepthScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_depth_scatter.pdf";
        String abdDepthStanScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_depth_stan_scatter.pdf";

        String dDepthDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_depth_density.pdf";
        String dDepthSDDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_SD_density.pdf";
        String dDepthStanDesity = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_depth_stan_density.pdf";
        String dDepthScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_depth_scatter.pdf";
        String dDepthStanScatter = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_depth_stan_scatter.pdf";

        this.plotPopDep(abSampleFileS, abDepthDesity, abDepthSDDesity, abDepthStanDesity, abDepthScatter, abDepthStanScatter, "AB");
        this.plotPopDep(abdSampleFileS, abdDepthDesity, abdDepthSDDesity, abdDepthStanDesity, abdDepthScatter, abdDepthStanScatter, "ABD");
        this.plotPopDep(dSampleFileS, dDepthDesity, dDepthSDDesity, dDepthStanDesity, dDepthScatter, dDepthStanScatter, "D");
    }

    private void plotPopDep (String infileS, String depthDensity, String sdDensity, String depthStanDensity, String depthScatter, String depthStanScatter, String genomeType) {
        RowTable<String> t = new RowTable<>(infileS);
        double[] depth = t.getColumnAsDoubleArray(1);
        double[] depthSD = t.getColumnAsDoubleArray(2);
        double[] depthStan = t.getColumnAsDoubleArray(3);
        double[] depthStanSD = t.getColumnAsDoubleArray(4);
        DensityPlot d = new DensityPlot(depth);
        d.setTitle(genomeType);
        d.setXLim(0, 20);
        d.setXLab("Mean of depth");
        d.setYLab("Density");
        d.saveGraph(depthDensity);

        d = new DensityPlot(depthSD);
        d.setTitle(genomeType);
        d.setXLim(0, 13);
        d.setXLab("Mean of depth SD");
        d.setYLab("Density");
        d.setSmoothN(5000);
        //d.showGraph();
        d.saveGraph(sdDensity);

        d = new DensityPlot(depthStan);
        d.setXLim(0,3);
        d.setXLab("Mean of standardized depth");
        d.setYLab("Density");
        d.saveGraph(depthStanDensity);
        ScatterPlot s = new ScatterPlot(depth, depthSD);
        s.setTitle(genomeType);
        if (genomeType.equals("AB")) {
            s.setXLim(0, 8);
            s.setYLim(0, 8);
        }
        else {
            s.setXLim(0, 20);
            s.setYLim(0, 12);
        }

        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.setColor(255, 0, 0, 20);
        s.saveGraph(depthScatter);
        s = new ScatterPlot(depthStan, depthStanSD);
        s.setXLim(0, 2);
        s.setYLim(0, 1);
        s.setXLab("Mean of standardized depth");
        s.setYLab("SD of standardized depth");
        s.setColor(255, 0, 0, 20);
        s.saveGraph(depthStanScatter);
    }

    public void popdepPlotSample() {
        String abFile1 = "/Volumes/VMap2_Fei/popdep_vmap2/AB/chr001_AB_popdep_vmap2.txt.gz";
        String abFile2 = "/Volumes/VMap2_Fei/popdep_vmap2/AB/chr002_AB_popdep_vmap2.txt.gz";
        String abdFile1 = "/Volumes/VMap2_Fei/popdep_vmap2/ABD/chr001_ABD_popdep_vmap2.txt.gz";
        String abdFile2 = "/Volumes/VMap2_Fei/popdep_vmap2/ABD/chr002_ABD_popdep_vmap2.txt.gz";
        String dFile1 = "/Volumes/VMap2_Fei/popdep_vmap2/D/chr005_D_popdep_vmap2.txt.gz";
        String dFile2 = "/Volumes/VMap2_Fei/popdep_vmap2/D/chr006_D_popdep_vmap2.txt.gz";
        int sampleSize = 50000;
        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_sample.txt";
        this.samplePopDep(abFile1, abFile2, abSampleFileS, sampleSize);
        this.samplePopDep(abdFile1, abdFile2, abdSampleFileS, sampleSize);
        this.samplePopDep(dFile1, dFile2, dSampleFileS, sampleSize);
    }

    private void samplePopDep (String file1, String file2, String outfileS, int sampleSize) {
        File[] inFiles = new File[2];
        inFiles[0] = new File(file1);
        inFiles[1] = new File(file2);
        int[] chrs = new int[inFiles.length];
        int[] chrLengths = new int[2];
        int[] sizes = new int[inFiles.length];
        int tLength = 0;
        for (int i = 0; i < chrs.length; i++) {
            chrs[i] = Integer.parseInt(inFiles[i].getName().split("_")[0].replaceFirst("chr", ""));
            chrLengths[i] = RefV1Utils.getChrIDLength(chrs[i]);
            tLength+=chrLengths[i];
        }
        String header = null;
        List<String> cList = new ArrayList<>();
        for (int i = 0; i < chrs.length; i++) {
            sizes[i] = (int)Math.round(sampleSize*((double)chrLengths[i]/tLength));
            int step = chrLengths[i]/sizes[i];
            try {
                int cnt = -1;
                BufferedReader br = IOUtils.getTextGzipReader(inFiles[i].getAbsolutePath());
                header = br.readLine();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%step != 0) {
                        if (cnt%10000000 == 0) {
                            System.out.println(cnt);
                        }
                        continue;
                    }
                    cList.add(temp);

                }
                System.out.println(inFiles[i].getName());
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            List<String> l = new ArrayList<>();
            for (int i = 0; i < cList.size(); i++) {
                bw.write(cList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {

        }


    }

    public void runStep2 () {
        String parameterDirS = "/data1/home/feilu/popdepStep2/parameters";
        List<File> fList = IOUtils.getFileListInDirEndsWith(parameterDirS, ".txt");
        for (int i = 0; i < fList.size(); i++) {
            new PopDep(fList.get(i).getAbsolutePath());
            System.out.println(String.valueOf(i+1)+" "+fList.get(i).getName());
        }
        System.out.println("Done");
    }

    public void mkStep2FileS() {
        String sourceFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/inputfile/parameters_popdep_2.txt";
        String dTaxaModeFileS = "/data1/home/feilu/popdep/dTaxaDepthMode.txt";
        String abTaxaModeFileS = "/data1/home/feilu/popdep/abTaxaDepthMode.txt";
        String abdTaxaModeFileS = "/data1/home/feilu/popdep/abdTaxaDepthMode.txt";
        String dTaxaBamFileS = "/data1/home/feilu/popdep/dTaxaBam.txt";
        String abTaxaBamFileS = "/data1/home/feilu/popdep/abTaxaBam.txt";
        String abdTaxaBamFileS = "/data1/home/feilu/popdep/abdTaxaBam.txt";
        String samPath = "/data1/programs/samtools-1.8/samtools";
        String dOutDirS = "/data1/home/feilu/popdepStep2/D";
        String abOutDirS = "/data1/home/feilu/popdepStep2/AB";
        String abdOutDirS = "/data1/home/feilu/popdepStep2/ABD";
        String step2DirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step2";
        List<String> template = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(sourceFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                template.add(temp);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.mkStep2(template, dTaxaBamFileS, dTaxaModeFileS, samPath, dOutDirS, step2DirS, "D");
        this.mkStep2(template, abTaxaBamFileS, abTaxaModeFileS, samPath, abOutDirS, step2DirS, "AB");
        this.mkStep2(template, abdTaxaBamFileS, abdTaxaModeFileS, samPath, abdOutDirS, step2DirS, "ABD");
    }

    private void mkStep2 (List<String> template, String taxaBamFileS, String taxaModeFileS, String samPath, String outDirS, String step2, String genomeType) {
        StringBuilder sb = new StringBuilder();
        TIntArrayList indexList = new TIntArrayList();
        for (int i = 0; i < template.size(); i++) {
            String current = template.get(i);
            if (current.startsWith("@")) continue;
            else if (current.startsWith("#")) continue;
            else if (current.equals("")) continue;
            indexList.add(i);
        }
        int[] indices = indexList.toArray(new int[indexList.size()]);
        int[] chrIDs = null;
        TIntArrayList chrIDList = new TIntArrayList();
        if (genomeType.equals("AB")) {
            chrIDList.addAll(RefV1Utils.getChrIDsOfSubgenomeA());
            chrIDList.addAll(RefV1Utils.getChrIDsOfSubgenomeB());
        }
        else if (genomeType.equals("ABD")) {
            chrIDList.addAll(RefV1Utils.getChrIDs());
        }
        else if (genomeType.equals("D")) {
            chrIDList.addAll(RefV1Utils.getChrIDsOfSubgenomeD());
        }
        chrIDs = chrIDList.toArray();
        Arrays.sort(chrIDs);
        for (int i = 0; i < chrIDs.length; i++) {
            sb.setLength(0);
            sb.append("chr").append(PStringUtils.getNDigitNumber(3, chrIDs[i])).append("_").append(genomeType).append("_popdep_vmap2.txt.gz");
            String outfileS = new File(outDirS,sb.toString()).getAbsolutePath();
            sb.setLength(0);
            sb.append(genomeType).append("_").append(chrIDs[i]).append("_parameters_popdep_2.txt");
            String step2FileS = new File(step2, sb.toString()).getAbsolutePath();
            try {
                int cnt = 0;
                BufferedWriter bw = IOUtils.getTextWriter(step2FileS);
                for (int j = 0; j < template.size(); j++) {
                    int index = Arrays.binarySearch(indices, j);
                    if (index < 0) {
                        bw.write(template.get(j));
                    }
                    else {
                        if (cnt == 0) {
                            bw.write(taxaBamFileS);
                        }
                        else if (cnt == 1) {
                            bw.write(taxaModeFileS);
                        }
                        else if (cnt == 2) {
                            bw.write(String.valueOf(chrIDs[i]));
                        }
                        else if (cnt == 3) {
                            bw.write(String.valueOf(RefV1Utils.getChrIDLength(chrIDs[i])));
                        }
                        else if (cnt == 4) {
                            bw.write("5");
                        }
                        else if (cnt == 5) {
                            bw.write(samPath);

                        }
                        else if (cnt == 6) {
                            bw.write("32");
                        }
                        else if (cnt == 7) {
                            bw.write(outfileS);
                        }
                        cnt++;
                    }
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

    public void runStep1 () {
        String abdInfileS = "/data1/home/feilu/popdep/abd_parameters_popdep_1.txt";
        String abInfileS = "/data1/home/feilu/popdep/ab_parameters_popdep_1.txt";
        String dInfileS = "/data1/home/feilu/popdep/d_parameters_popdep_1.txt";
        //new PopDep(dInfileS);
//        new PopDep(abInfileS);
//        new PopDep(abdInfileS);
    }

    public void mkChrLengthFile () {
        String abdLengthFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/abdChromosomeLength.txt";

        try {
            BufferedWriter bw = IOUtils.getTextWriter(abdLengthFileS);
            bw.write("Chromosome(Should be all choromosomes)\tLength");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            int[] ids = RefV1Utils.getChrIDs();
            for (int i = 0; i < RefV1Utils.getChrIDs().length; i++) {
                sb.setLength(0);
                sb.append(ids[i]).append("\t").append(RefV1Utils.getChrIDLength(ids[i]));
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

    public void mkTaxaBamFile () {
        String vmap2TaxaFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/wheatVMapII_germplasmInfo_20191225.txt";
        String old2newFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/bam_Old2New.map.txt";
        String abdFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/002_taxaRefBam.ABDgenome.manual.addNAFU.txt";
        String abFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/003_taxaRefBam.Dgenome.addNAFU.txt";
        String dFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/004_taxaRefBam.ABgenome.removeBadTaxa.addNAFU.addS1.txt";
        String abdOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/abdTaxaBam.txt";
        String abOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/abTaxaBam.txt";
        String dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/dTaxaBam.txt";
        RowTable<String> t = new RowTable<>(old2newFileS);
        HashMap<String, String> old2new = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCell(i,0).endsWith(".bai")) continue;
            old2new.put(t.getCell(i,0), t.getCell(i,1));
        }
        t = new RowTable<>(vmap2TaxaFileS);
        HashSet<String> set = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.getCell(i,25).equals("1")) continue;
            set.add(t.getCell(i, 4));
        }
        String[] taxa = set.toArray(new String[set.size()]);
        Arrays.sort(taxa);
        System.out.println(taxa.length);
        this.mkTaxaBam(old2new, taxa, abdFileS, abdOutfileS);
        this.mkTaxaBam(old2new, taxa, abFileS, abOutfileS);
        this.mkTaxaBam(old2new, taxa, dFileS, dOutfileS);

    }

    private void mkTaxaBam (HashMap<String, String> old2new, String[] taxa, String infileS, String outfileS) {
        RowTable<String> t = new RowTable<>(infileS);
        List<String> taxaList = t.getColumn(0);
        HashSet<String> set = new HashSet<>(taxaList);
        set.retainAll(Arrays.asList(taxa));
        String[] nTaxa = set.toArray(new String[set.size()]);
        Arrays.sort(nTaxa);
        List<String>[] bamLists = new ArrayList[nTaxa.length];
        for (int i = 0; i < nTaxa.length; i++) {
            bamLists[i] = new ArrayList<>();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(nTaxa, t.getCell(i, 0));
            if (index < 0) continue;
            bamLists[index].add(old2new.get(t.getCell(i,2)));
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa	Bams(A list of bams of the taxon, seperated by the delimiter of Tab)");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < nTaxa.length; i++) {
                sb.setLength(0);
                sb.append(nTaxa[i]);
                for (int j = 0; j < bamLists[i].size(); j++) {
                    sb.append("\t").append(bamLists[i].get(j));
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
}
