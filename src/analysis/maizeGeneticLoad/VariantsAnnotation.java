/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maizeGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author feilu
 */
class VariantsAnnotation {
    
    public VariantsAnnotation () {
        this.mkHeaderOfHap3();
    }
    
    public void mkHeaderOfHap3 () {
        String infileDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/hmp321_agp4";
        String outfileDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/hmp/hmp321_header";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String oufileS = new File(outfileDirS, f.getName().replaceFirst(".gz", ".header.txt")).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(oufileS);
                bw.write("Chr\tPos\tRef\tAlt\tMaf");
                bw.newLine();
                String temp = null;
                List<String> l = null;
                StringBuilder sb = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    cnt++;
                    if (cnt%100000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    sb = new StringBuilder();
                    temp = temp.substring(0, 160);
                    l = PStringUtils.fastSplit(temp);
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    l = PStringUtils.fastSplit(l.get(7), "MAF=");
                    l = PStringUtils.fastSplit(l.get(1), ";");
                    sb.append(l.get(0));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
}
