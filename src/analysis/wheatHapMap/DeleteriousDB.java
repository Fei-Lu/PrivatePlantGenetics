/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import format.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import graphcis.r.DensityPlot;
import graphcis.r.ScatterPlot;
import graphcis.tablesaw.TablesawUtils;
import static org.apache.commons.math3.stat.inference.TestUtils.t;
import tech.tablesaw.api.Table;
import tech.tablesaw.io.csv.CsvReadOptions;

/**
 *
 * @author feilu
 */
public class DeleteriousDB {
    
    public DeleteriousDB () {
        this.selectHCGenes();
    }
    
    public void selectHCGenes () {
        String infileS = "/Users/feilu/Documents/database/wheat/gene/gene_expression/geneExpression.txt";
        double tpmThresh = 0.1;
        
        try {
            Table t = TablesawUtils.readTsv(infileS);
            int a = 3;
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        
//        TDoubleArrayList meanList = new TDoubleArrayList();
//        TDoubleArrayList sdList = new TDoubleArrayList();
//        RowTable<String> t = new RowTable(infileS);
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            double mean = t.getCellAsDouble(i, 5);
//            if (mean < tpmThresh) continue;
//            meanList.add(mean);
//            sdList.add(t.getCellAsDouble(i, 6));
//        }
//        System.out.println(meanList.size());
//        double[] means = meanList.toArray();
//        double[] sds = sdList.toArray();
//        
//        ScatterPlot s = new ScatterPlot(means, sds);
//        //s.showGraph();
//        DensityPlot d = new DensityPlot(means);
//        d.showGraph();
        
    }
}
