/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.VMap2;

import pgl.infra.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.graph.r.Histogram;
import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
public class AnnoSift {
    
    public AnnoSift() {
        this.genicStats();
    }
    
    public void genicStats () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation";
        String outpdf = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/annoGenicStats/siftDistribution.pdf";
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        TDoubleArrayList siftList = new TDoubleArrayList();
        HashSet<String> typeSet = new HashSet();
        double siftThresh = 0.05;
        int sampleSize = 10000;
        fList.stream().forEach(f -> {
            RowTable<String> t = new RowTable (f.getAbsolutePath());
            int cIndex = t.getColumnIndex("SIFT_score");
            String cell = null;
            String type = null;
            double value = -1;
            for (int i = 0; i < t.getRowNumber(); i++) {
                cell = t.getCell(i, cIndex);
                type = t.getCell(i, cIndex-1);
                if (cell.startsWith("N")) continue;
                value = Double.parseDouble(cell);
                siftList.add(value);
            }
            
        });
        System.out.println(siftList.size());
        int cnt = 0;
        for (int i = 0; i < siftList.size(); i++) {
            if (siftList.get(i) > 0.05) continue;
            cnt++;
        }
        System.out.println(cnt);
        siftList.shuffle(new Random());
        
        double[] siftValues = siftList.toArray(0, sampleSize);
        
        Histogram h = new Histogram(siftValues);
        h.setTitle("Sift value distribution of VMap II");
        h.setXLab("Sift value");
        h.setYLab("Proportion");
        h.saveGraph(outpdf);
    }
    
}
