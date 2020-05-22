/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.data;

import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
public class WheatReferenceGene {
    
    public WheatReferenceGene () {
        //this.GFFToPGF();
        //this.mkNonoverlappingGene();
        //this.mkExpressionByTissue();
        //this.mkExpressionSummary();
//       this.mkExpressedGene();

    }
    
    public void mkExpressedGene () {
        String pgfFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String infileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/expressionByTissueSummary.txt";
        String outfileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/geneExpression.txt";
        String header = "Gene\tIs_Overlapped_gene(1,0)\tNumber_Overlapped_gene\tIs_Unique_gene(1,0)\tLongest_transcript\tTPM_mean\tTPM_sd\tTPM_rsd\tTPM_max\tTPM_min";
        RowTable<String> t = new RowTable(infileS);
        t.sortAsText(0);
        List<String> tList = t.getColumn(0);
        GeneFeature gf = new GeneFeature (pgfFileS);
        Range[] geneRanges = new Range[gf.getGeneNumber()];
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            geneRanges[i] = new Range(gf.getGeneChromosome(i), gf.getGeneStart(i), gf.getGeneEnd(i));
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            TIntArrayList indexList = new TIntArrayList();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < geneRanges.length; i++) {
                indexList.clear();
                sb.setLength(0);
                sb.append(gf.getGeneName(i)).append("\t");
                for (int j = 0; j < geneRanges.length; j++) {
                    if (i == j) continue;
                    if (geneRanges[i].getRangeChromosome() != geneRanges[j].getRangeChromosome()) continue;
                    if (geneRanges[i].isOverlap(geneRanges[j])) {
                        indexList.add(j);
                    }
                }
                if (indexList.size() == 0) {
                    sb.append(0).append("\t").append(0).append("\t").append(1).append("\t");
                }
                else {
                    sb.append(1).append("\t").append(indexList.size()).append("\t");
                    int value = 1;
                    for (int j = 0; j < indexList.size(); j++) {
                        if (geneRanges[i].getRangeSize() < geneRanges[indexList.get(j)].getRangeSize()) {
                            value = 0;
                            break;
                        }
                    }
                    sb.append(value).append("\t");
                }
                String tName = gf.getTranscriptName(i, gf.getLongestTranscriptIndex(i));
                sb.append(tName).append("\t");
                int index = Collections.binarySearch(tList, tName);
                sb.append(t.getCell(index, 1)).append("\t").append(t.getCell(index, 2)).append("\t").append(t.getCell(index, 3)).append("\t");
                sb.append(t.getCell(index, 4)).append("\t").append(t.getCell(index, 5));
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
    
    
    public void mkExpressionSummary () {
        String infileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/expressionByTissue.txt";
        String outfileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/expressionByTissueSummary.txt";
        RowTable<String> t = new RowTable(infileS);
        String header = "Transcripts\tTPM_mean\tTPM_sd\tTPM_rsd\tTPM_max\tTPM_min";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            double[] values = new double[t.getColumnNumber()-1];
            StringBuilder sb = new StringBuilder();
            DescriptiveStatistics ds = null;
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append(t.getCell(i, 0));
                for (int j = 0; j < values.length; j++) {
                    values[j] = t.getCellAsDouble(i, j+1);
                }
                ds = new DescriptiveStatistics(values);
                double mean = ds.getMean();
                double sd = ds.getStandardDeviation();
                double rsd = sd/mean;
                double max = ds.getMax();
                double min = ds.getMin();
                sb.append("\t").append((float)mean).append("\t").append((float)sd).append("\t").append((float)rsd);
                sb.append("\t").append((float)max).append("\t").append(min);
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
    
    public void mkExpressionByTissue () {
        String sampleInFileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/download/Table_S1.txt";
        String expressionInFileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/download/Development_tpm.tsv.gz";
        String outfileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/expressionByTissue.txt";
        RowTable<String> t = new RowTable(sampleInFileS);
        HashSet<String> s = new HashSet(t.getColumn(9));
        String[] tissues = s.toArray(new String[s.size()]);
        Arrays.sort(tissues);
        List<String>[] sampleLists = new List[tissues.length];
        for (int i = 0; i < sampleLists.length; i++) sampleLists[i] = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(tissues, t.getCell(i, 9));
            sampleLists[index].add(t.getCell(i, 0));
        }
        int[][] sampleIndex = new int[sampleLists.length][];
        for (int i = 0; i < sampleLists.length; i++) {
            sampleIndex[i] = new int[sampleLists[i].size()];
            Collections.sort(sampleLists[i]);
        }
        t = new RowTable(expressionInFileS);
        List<String> header = t.getHeader();
        for (int i = 0; i < header.size(); i++) {
            for (int j = 0; j < tissues.length; j++)  {
                int index = Collections.binarySearch(sampleLists[j], header.get(i));
                if (index < 0) continue;
                sampleIndex[j][index] = i;
            }
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("Transcripts");
            for (int i = 0; i < tissues.length; i++) {
                sb.append("\t").append(tissues[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append(t.getCell(i, 0));
                for (int j = 0; j < tissues.length; j++) {
                    double e = 0;
                    for (int k = 0; k < sampleLists[j].size(); k++) {
                        e+=Double.parseDouble(t.getCell(i, sampleIndex[j][k]));
                    }
                    e = e/sampleIndex[j].length;
                    sb.append("\t").append((float)e);
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

    public void mkNonoverlappingGff3 () {
        String infileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_nonoverlap.txt";
        String gff3FileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3";
        String outGff3FileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab_nonoverlap.gff3";
        RowTable<String> t = new RowTable<> (infileS);
        
    }

    public void mkNonoverlappingGene () {
        String pgfFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String outfileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_nonoverlap.txt";
        String header = "Gene\tIs_Overlapped_gene(1,0)\tNumber_Overlapped_gene\tIs_Unique_gene(1,0)\tLongest_transcript";
        GeneFeature gf = new GeneFeature (pgfFileS);
        Range[] geneRanges = new Range[gf.getGeneNumber()];
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            geneRanges[i] = new Range(gf.getGeneChromosome(i), gf.getGeneStart(i), gf.getGeneEnd(i));
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            TIntArrayList indexList = new TIntArrayList();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < geneRanges.length; i++) {
                indexList.clear();
                sb.setLength(0);
                sb.append(gf.getGeneName(i)).append("\t");
                for (int j = 0; j < geneRanges.length; j++) {
                    if (i == j) continue;
                    if (geneRanges[i].getRangeChromosome() != geneRanges[j].getRangeChromosome()) continue;
                    if (geneRanges[i].isOverlap(geneRanges[j])) {
                        indexList.add(j);
                    }
                }
                if (indexList.size() == 0) {
                    sb.append(0).append("\t").append(0).append("\t").append(1).append("\t");
                }
                else {
                    sb.append(1).append("\t").append(indexList.size()).append("\t");
                    int value = 1;
                    for (int j = 0; j < indexList.size(); j++) {
                        if (geneRanges[i].getRangeSize() < geneRanges[indexList.get(j)].getRangeSize()) {
                            value = 0;
                            break;
                        }
                    }
                    sb.append(value).append("\t");
                }
                String tName = gf.getTranscriptName(i, gf.getLongestTranscriptIndex(i));
                sb.append(tName);
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
    
    public void GFFToPGF () {
        String inputGFF = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3";
        String outputPGF = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature();
        gf.readFromWheatGFF(inputGFF);
        gf.writeFile(outputPGF);
    }
    
}
