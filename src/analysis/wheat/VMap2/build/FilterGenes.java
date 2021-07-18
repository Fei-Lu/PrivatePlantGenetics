package analysis.wheat.VMap2.build;

import com.koloboke.collect.map.hash.HashIntFloatMap;
import com.koloboke.collect.map.hash.HashIntFloatMaps;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.graph.r.ScatterPlot;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

class FilterGenes {
    public FilterGenes() {
        this.geneDepthProfile();
//        this.geneDepthPlot();
//        this.geneDepthPCAPlot();
//        this.densityFilter();
    }

    public void densityFilter() {
        String infileS = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/GenePopDepth_calByMergedSub.txt";
        String uniqueFileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/geneExpression.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/syntenic_genes.txt";

        RowTable<String> t = new RowTable<>(uniqueFileS);
        HashMap<String, String> uniqueMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            uniqueMap.put(t.getCell(i,0), t.getCell(i, 3));
        }
        t = new RowTable<>(infileS);
        DoubleArrayList[] abList = new DoubleArrayList[2];
        DoubleArrayList[] dList = new DoubleArrayList[2];
        for (int i = 0; i < 2; i++) {
            abList[i] = new DoubleArrayList();
            dList[i] = new DoubleArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            String type = t.getCell(i, 8);
            if (type.startsWith("D")) {
                dList[0].add(Double.parseDouble(t.getCell(i,10)));
                dList[1].add(Double.parseDouble(t.getCell(i,11)));
            }
            else {
                abList[0].add(Double.parseDouble(t.getCell(i,10)));
                abList[1].add(Double.parseDouble(t.getCell(i,11)));
            }
        }
        int binNum = 50;
        double proportion = 0.8;
        Grid grAB = new Grid(-0.3, 0.3, -0.3, 0.3, binNum);
        for (int i = 0; i < abList[0].size(); i++) {
            grAB.addXY(abList[0].getDouble(i), abList[1].getDouble(i));
        }
        grAB.buildHashMap();
        int abIndexThresh = grAB.getOrderIndexOfProportionOfSite(proportion);
        Grid grD = new Grid(-0.15, -0.0, -0.08, 0.08, binNum);
        for (int i = 0; i < dList[0].size(); i++) {
            grD.addXY(dList[0].getDouble(i), dList[1].getDouble(i));
        }
        grD.buildHashMap();
        int dIndexThresh = grD.getOrderIndexOfProportionOfSite(proportion);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("GeneName	Chr	PosStart	PosEnd\tIs_Unique_gene(1,0)");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                String type = t.getCell(i, 8);
                double x = Double.parseDouble(t.getCell(i,10));
                double y = Double.parseDouble(t.getCell(i,11));
                boolean pass = false;
                if (type.startsWith("D")) {
                    if (grD.isHighDensity(x, y , dIndexThresh)) {
                        pass = true;
                    }
                }
                else {
                    if (grAB.isHighDensity(x, y , abIndexThresh)) {
                        pass = true;
                    }
                }
                if (pass == false) continue;
                sb.setLength(0);
                sb.append(t.getCell(i,0)).append("\t").append(t.getCell(i,1)).append("\t").append(t.getCell(i,2))
                        .append("\t").append(t.getCell(i,3)).append("\t").append(uniqueMap.get(t.getCell(i,0)));
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

    public void geneDepthPCAPlot() {
        //Aoyue did the PCA, for AB, PC1 = 92.83%, PC2 = 5.05%; for D, PC1 = 87.94%, PC2 = 11.59%
        String infileS = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/GenePopDepth_calByMergedSub.txt";
        String outfileSAB = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/GeneDepthFeature_AB.pdf";
        String outfileSD = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/GeneDepthFeature_D.pdf";
        String fOutfileSAB = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/GeneDepthFeature_AB_filter.pdf";
        String fOutfileSD = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/GeneDepthFeature_D_filter.pdf";
        RowTable<String> t = new RowTable<>(infileS);
        DoubleArrayList[] abList = new DoubleArrayList[2];
        DoubleArrayList[] dList = new DoubleArrayList[2];
        for (int i = 0; i < 2; i++) {
            abList[i] = new DoubleArrayList();
            dList[i] = new DoubleArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            String type = t.getCell(i, 8);
            if (type.startsWith("D")) {
                dList[0].add(Double.parseDouble(t.getCell(i,10)));
                dList[1].add(Double.parseDouble(t.getCell(i,11)));
            }
            else {
                abList[0].add(Double.parseDouble(t.getCell(i,10)));
                abList[1].add(Double.parseDouble(t.getCell(i,11)));
            }
        }
        ScatterPlot s = new ScatterPlot(abList[0].toDoubleArray(), abList[1].toDoubleArray());
        s.setTitle("GeneDepthFeature_AB");
        s.setColor(255, 0, 0, 10);
        s.setXLim(-0.5,0.5);
        s.setYLim(-0.5,0.5);
        s.setXLab("PC1:92.83%");
        s.setYLab("PC2:5.05%");
        s.saveGraph(outfileSAB);
        s = new ScatterPlot(dList[0].toDoubleArray(), dList[1].toDoubleArray());
        s.setTitle("GeneDepthFeature_D");
        s.setColor(255, 0, 0, 10);
        s.setXLim(-0.3,0.3);
        s.setYLim(-0.3,0.3);
        s.setXLab("PC1:87.94%");
        s.setYLab("PC2:11.59%");
        s.saveGraph(outfileSD);

        int binNum = 50;
        double proportion = 0.8;
        Grid gr = new Grid(-0.3, 0.3, -0.3, 0.3, binNum);
        for (int i = 0; i < abList[0].size(); i++) {
            gr.addXY(abList[0].getDouble(i), abList[1].getDouble(i));
        }
        gr.buildHashMap();
        DoubleArrayList xList = new DoubleArrayList();
        DoubleArrayList yList = new DoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportion);
        int cnt = 0;
        for (int i = 0; i < abList[0].size(); i++) {
            double x = abList[0].getDouble(i);
            double y = abList[1].getDouble(i);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
            cnt++;
        }
        System.out.println((double)cnt/abList[0].size());
        s = new ScatterPlot(xList.toDoubleArray(), yList.toDoubleArray());
        s.setTitle("GeneDepthFeature_filter_AB");
        s.setColor(255, 0, 0, 10);
        s.setXLim(-0.5,0.5);
        s.setYLim(-0.5,0.5);
        s.setXLab("PC1:92.83%");
        s.setYLab("PC2:5.05%");
        s.saveGraph(fOutfileSAB);

        gr = new Grid(-0.15, -0.0, -0.08, 0.08, binNum);
        for (int i = 0; i < dList[0].size(); i++) {
            gr.addXY(dList[0].getDouble(i), dList[1].getDouble(i));
        }
        gr.buildHashMap();
        xList = new DoubleArrayList();
        yList = new DoubleArrayList();
        indexThresh = gr.getOrderIndexOfProportionOfSite(proportion);
        cnt = 0;
        for (int i = 0; i < dList[0].size(); i++) {
            double x = dList[0].getDouble(i);
            double y = dList[1].getDouble(i);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
            cnt++;
        }
        System.out.println((double)cnt/dList[0].size());
        s = new ScatterPlot(xList.toDoubleArray(), yList.toDoubleArray());
        s.setTitle("GeneDepthFeature_filter_D");
        s.setColor(255, 0, 0, 10);
        s.setXLim(-0.3,0.3);
        s.setYLim(-0.3,0.3);
        s.setXLab("PC1:87.94%");
        s.setYLab("PC2:11.59%");
        s.saveGraph(fOutfileSD);
    }

    public void geneDepthPlot () {
        String infileS = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/geneDepth.txt";
        DoubleArrayList[] abLists = new DoubleArrayList[4];
        DoubleArrayList[] dLists = new DoubleArrayList[4];
        for (int i = 0; i < abLists.length; i++) {
            abLists[i] = new DoubleArrayList();
            dLists[i] = new DoubleArrayList();
        }
        int[] aChrIDs = RefV1Utils.getChrIDsOfSubgenomeA();
        int[] bChrIDs = RefV1Utils.getChrIDsOfSubgenomeB();
        int[] dChrIDs = RefV1Utils.getChrIDsOfSubgenomeD();
        int[] abChrIDs = new int[aChrIDs.length+bChrIDs.length];
        System.arraycopy(aChrIDs,0, abChrIDs,0, aChrIDs.length);
        System.arraycopy(bChrIDs,0, abChrIDs, aChrIDs.length, bChrIDs.length);
        RowTable<String> t = new RowTable<>(infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chr = Integer.parseInt(t.getCell(i,1));
            if (chr == 0) continue;
            if (Arrays.binarySearch(abChrIDs, chr) >=0) {
                for (int j = 0; j < 4; j++) {
                    abLists[j].add(Double.parseDouble(t.getCell(i, j+4)));
                }
            }
            else if (Arrays.binarySearch(dChrIDs,chr) >= 0) {
                for (int j = 0; j < 4; j++) {
                    dLists[j].add(Double.parseDouble(t.getCell(i, j+4)));
                }
            }

        }

        ScatterPlot s = new ScatterPlot(abLists[0].toDoubleArray(), abLists[2].toDoubleArray());
        s.setColor(255, 0, 0, 2);
        s.setXLim(0,20);
        s.setYLim(0,20);
        s.showGraph();

    }

    public void geneDepthProfile () {
        String geneFeatureFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String abPopDepDirS = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/AB";
        String dPopDepDirS = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/D";
        String outfileS = "/Users/feilu/Documents/analysisL/production/vmap2/filterGenes/geneDepth.txt";
        List<File> fList = IOUtils.getFileListInDirEndsWith(abPopDepDirS, ".gz");
        List<File> dfList = IOUtils.getFileListInDirEndsWith(dPopDepDirS, ".gz");
        fList.addAll(dfList);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        HashIntFloatMap posDepthMap = HashIntFloatMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        HashIntFloatMap posDepthSDMap = HashIntFloatMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("GeneName\tChr\tPosStart\tPosEnd\tSiteDepth\tSiteDepthSDBySites\tSiteDepthSD\tSiteDepthSDSDBySites");
            bw.newLine();
            int currentChr = -1;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < gf.getGeneNumber(); i++) {
                sb.setLength(0);
                sb.append(gf.getGeneName(i)).append("\t").append(gf.getChromosomeOfGene(i)).append("\t").append(gf.getGeneStart(i)).append("\t").append(gf.getGeneEnd(i)).append("\t");
                int chr = gf.getChromosomeOfGene(i);
                if (currentChr != chr) {
                    File currentFile = this.getPopDepPath(chr, fList);
                    if (currentFile == null) {
                        sb.append("NA\tNA\tNA\tNA");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    else {
                        posDepthMap.clear();
                        posDepthSDMap.clear();
                        BufferedReader br = IOUtils.getTextGzipReader(currentFile.getAbsolutePath());
                        String temp = br.readLine();
                        List<String> l = null;
                        while ((temp = br.readLine()) != null) {
                            l = PStringUtils.fastSplit(temp);
                            posDepthMap.put(Integer.parseInt(l.get(0)), Float.parseFloat(l.get(1)));
                            posDepthSDMap.put(Integer.parseInt(l.get(0)), Float.parseFloat(l.get(2)));
                        }
                        br.close();
                        System.out.println("Finished reading in " + currentFile.getName());
                        currentChr = chr;
                    }
                }
                int tIndex = gf.getLongestTranscriptIndex(i);
                List<Range> cdsList = gf.getCDSList(i, tIndex);
                DoubleArrayList depthList = new DoubleArrayList();
                DoubleArrayList depthSDList = new DoubleArrayList();
                for (int j = 0; j < cdsList.size(); j++) {
                    int cdsLength = cdsList.get(j).end - cdsList.get(j).start;
                    for (int k = 0; k < cdsLength; k++) {
                        depthList.add(posDepthMap.get(cdsList.get(j).start+k));
                        depthSDList.add(posDepthSDMap.get(cdsList.get(j).start+k));
                    }
                }
                DescriptiveStatistics d = new DescriptiveStatistics(depthList.toDoubleArray());
                sb.append((float)d.getMean()).append("\t").append((float)d.getStandardDeviation()).append("\t");
                d = new DescriptiveStatistics(depthSDList.toDoubleArray());
                sb.append((float)d.getMean()).append("\t").append((float)d.getStandardDeviation());
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

    private File getPopDepPath(int chr, List<File> fList) {
        String header = "chr"+PStringUtils.getNDigitNumber(3, chr);
        for (int i = 0; i < fList.size(); i++) {
            if (fList.get(i).getName().startsWith(header)) return fList.get(i);
        }
        return null;
    }
}
