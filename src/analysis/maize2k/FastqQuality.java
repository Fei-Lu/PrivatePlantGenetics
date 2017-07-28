/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maize2k;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashSet;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class FastqQuality {
    
    public FastqQuality () {
        this.sampleReads();
        
    }
    
    private void sampleReads () {
        String infileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/rawdata";
        String outputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/maize2k/sampleSeq";
        int readNum = 10000;
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        int a = 3;
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfile1 = new File (outputDirS, name+"_1.fq.gz").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.fq.gz").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                String temp = null;
                int n = 0;
                while ((temp = br1.readLine()) != null) {
                    n++;
                    br1.readLine();br1.readLine();br1.readLine();
                    if (n%1000000 == 0) {
                        //System.out.println(String.valueOf(n)+"\t"+infile1);
                    }
                }
                br1.close();
                System.out.println(String.valueOf(n)+"\t"+infile1);
                int interval = n/readNum;
                
                br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfile2);
                String temp1 = null;
                String temp2 = null;
                for (int i = 0; i < n; i++) {
                    if (i%interval == 0) {
                        bw1.write(br1.readLine());bw1.write(br1.readLine());bw1.write(br1.readLine());bw1.write(br1.readLine());
                        bw2.write(br2.readLine());bw2.write(br2.readLine());bw2.write(br2.readLine());bw2.write(br2.readLine());
                    }
                    else {
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                    }
                }
                bw1.flush();bw1.close();
                bw2.flush();bw2.close();
                br1.close();
                br2.close();
                System.out.println(name+ " completed");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            
        });
    }
}
