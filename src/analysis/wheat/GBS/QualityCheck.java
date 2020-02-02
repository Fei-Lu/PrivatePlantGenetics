/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.GBS;

import pgl.analysis.pipeline.grt.TagAnnotations;
import pgl.graphcis.r.DensityPlot;
import java.io.BufferedWriter;
import java.io.File;
import pgl.utils.IOUtils;

/**
 *
 * @author feilu
 */
class QualityCheck {
    
    public QualityCheck () {
        this.checkReadNumOfSample();
    }
    
    public void checkReadNumOfSample () {
        String dirS = "/Users/feilu/Documents/analysisL/production/wheatGBS/pipeOutput/tagsBySample";
        String outfileS = "/Users/feilu/Documents/analysisL/production/wheatGBS/qualityCheck/readNumPerSample.txt";
        String outfilePdf = "/Users/feilu/Documents/analysisL/production/wheatGBS/qualityCheck/readNumPerSample.pdf";
//        String dirS = "/Users/feilu/Documents/analysisL/production/qianhuaGBS/pipeOutput/tagsBySample";
//        String outfileS = "/Users/feilu/Documents/analysisL/production/qianhuaGBS/qualityCheck/readNumPerSample.txt";
//        String outfilePdf = "/Users/feilu/Documents/analysisL/production/qianhuaGBS/qualityCheck/readNumPerSample.pdf";
        File[] fs = new File(dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sbin");
        double[] nums = new double[fs.length];
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Sample\tReadNum");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                TagAnnotations tc = new TagAnnotations(fs[i].getAbsolutePath());
                nums[i] = tc.getTotalReadNumber();
                bw.write(fs[i].getName()+"\t"+String.valueOf(nums[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        DensityPlot d = new DensityPlot(nums);
        d.setXLab("Read num");
        d.setYLab("Density");
        d.setTitle("Distribution of read num per taxon");
        d.setXLim(0, 200000);
        d.saveGraph(outfilePdf);
    }
}
