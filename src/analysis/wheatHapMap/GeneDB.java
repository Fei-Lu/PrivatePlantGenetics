/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import format.genomeAnnotation.GeneFeature;
import graphcis.tablesaw.TablesawUtils;
import java.io.File;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Table;
import tech.tablesaw.plotly.Plot;
import tech.tablesaw.plotly.api.Histogram;
import tech.tablesaw.plotly.api.ScatterPlot;
import tech.tablesaw.plotly.components.Figure;

/**
 *
 * @author feilu
 */
public class GeneDB {
    
    public GeneDB () {
        this.selectHCGenes1();
        //this.selectHCGenes2();
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
