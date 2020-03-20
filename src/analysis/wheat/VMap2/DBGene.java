/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.VMap2;

import pgl.infra.genomeAnnotation.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.ColumnTable;
import pgl.graphcis.tablesaw.TablesawUtils;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Table;
import tech.tablesaw.plotly.Plot;
import tech.tablesaw.plotly.api.Histogram;
import tech.tablesaw.plotly.api.ScatterPlot;
import tech.tablesaw.plotly.components.Figure;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author feilu
 */
public class DBGene {
    
    public DBGene() {
//        this.selectHCGenes1();
//        this.addGFFInfo();
//        this.addSiftInfo();
    }

    //not done
    public void addSiftInfo () {
        String dbFileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation";
        Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(dbFileS);
        String header = two.getFirstElement();
        List<String> recordList = two.getSecondElement();
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        HashMap<String, String>[] infoMap = new HashMap[fList.size()];
        AtomicInteger aCnt = new AtomicInteger();
        fList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().split("_")[0].replaceFirst("chr", ""))-1;
            infoMap[chrIndex] = new HashMap<>();
            ColumnTable<String> t = new ColumnTable<>(f.getAbsolutePath());
            int columnIndex = t.getColumnIndex("Transcript");
            HashSet<String> tSet = new HashSet<>(t.getColumn(columnIndex));
            aCnt.addAndGet(tSet.size());

        });
        System.out.println(aCnt.intValue());
    }

    public void addGFFInfo () {
        String infileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        String geneFeatureFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        gf.sortGeneByName();
        Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(infileS);
        String header = two.getFirstElement();
        List<String> recordList = two.getSecondElement();
        String query = null;
        List<String> l = null;
        int index = -1;
        int tn = 0;
        for (int i = 0; i < recordList.size(); i++) {
            l = PStringUtils.fastSplit(recordList.get(i));
            query = l.get(0);
            index = gf.getGeneIndex(query);
            tn = gf.getTranscriptNumber(index);
            StringBuilder sb =  new StringBuilder();
            for (int j = 0; j < tn; j++) {
                if (l.get(1).equals(gf.getTranscriptName(index, j))) {
                    int len = 0;
                    sb.setLength(0);
                    List<Range> cds = gf.getCDSList(index, j);
                    for (int k = 0; k < cds.size(); k++) {
                        len+=cds.get(k).getRangeSize();
                    }
                    sb.append(recordList.get(i)).append("\t").append(cds.size()).append("\t").append(len);
                    recordList.set(i, sb.toString());
                }
            }
        }
        VMapDBUtils.writeDB(header+"\tCDSExonNumber\tCDSLength", recordList, infileS);
    }
    
    public void selectHCGenes1 () {
        String infileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/geneExpression.txt";
        String geneFeatureFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        String outfile_mean_sd = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC_mean_sd.html";
        String outfile_mean = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC_mean.html";
        double tpmThresh = 0.1;
        try {
            GeneFeature gf = new GeneFeature (geneFeatureFileS);
            Table t = TablesawUtils.readTsv(infileS);
            System.out.println(t.structure());
            System.out.println(t.rowCount());
            t = t.where(t.intColumn("Is_Unique_gene(1,0)").isEqualTo(1).and(t.doubleColumn(8).isGreaterThan(tpmThresh)));
            t.removeColumns("Is_Overlapped_gene(1,0)", "Number_Overlapped_gene", "Is_Unique_gene(1,0)");
            gf.sortGeneByName();
            int[] chr = new int[t.rowCount()];
            int[] start = new int[t.rowCount()];
            int[] end = new int[t.rowCount()];
            int[] strand = new int[t.rowCount()];
            String query = null;
            String tName = null;
            for (int i = 0; i < t.rowCount(); i++) {
                query = t.getString(i, 0);
                tName = t.getString(i, 1);
                int geneIndex = gf.getGeneIndex(query);
                for (int j = 0; j < gf.getTranscriptNumber(geneIndex); j++) {
                    if (!(tName.equals(gf.getTranscriptName(geneIndex, j)))) continue;
                    chr[i] = gf.getTranscriptChromosome(geneIndex, j);
                    start[i] = gf.getTranscriptStart(geneIndex, j);
                    end[i] = gf.getTranscriptEnd(geneIndex, j);
                    strand[i] = gf.getTranscriptStrand(geneIndex, j);
                }
            }
            IntColumn chrs = IntColumn.create("Chr", chr);
            IntColumn starts = IntColumn.create("TranStart", start);
            IntColumn ends = IntColumn.create("TranEnd", end);
            IntColumn strands = IntColumn.create("TranStrand", strand);
            t.insertColumn(2, chrs);
            t.insertColumn(3, starts);
            t.insertColumn(4, ends);
            t.insertColumn(5, strands);
            t = t.where(t.intColumn("Chr").isGreaterThan(0));
            TablesawUtils.writeTsv(t, outfileS);
            Figure f = ScatterPlot.create("HC gene expression", t, "TPM_mean", "TPM_sd");
            Plot.show(f, new File(outfile_mean_sd));
            f = Histogram.create("Mean of HC gene expression (log10)", t.numberColumn("TPM_mean").log10());
            Plot.show(f, new File(outfile_mean));
            System.out.println(t.shape());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
