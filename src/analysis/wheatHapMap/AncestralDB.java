/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import gnu.trove.set.TCharSet;
import gnu.trove.set.hash.TCharHashSet;
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
public class AncestralDB {
    
    public AncestralDB () {
        this.splitMaf();
    }
    
    public void splitMaf () {
        String inDirS = "/Volumes/Fei_HDD_Mac/Gerp/asAlle/";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/003_ancestral/bySubgenome/";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".maf");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = f.getName().replaceFirst(".maf", ".aln");
            outfileS = new File(outDirS, outfileS).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String temp = null;
                String header = "Chr\tPos\tRef\tUrartu\tBarley";
                bw.write(header);
                bw.newLine();
                br.readLine();
                String seq1 = null;
                String seq2 = null;
                String seq3 = null;
                int pos = 0;
                List<String> l = null;
                String chr = null;
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    seq1 = br.readLine();seq2 = br.readLine();seq3 = br.readLine();br.readLine();
                    l = PStringUtils.fastSplitOnWhitespace(seq1);
                    chr = l.get(1).split("\\.")[1];
                    pos = Integer.parseInt(l.get(2))+1;
                    seq1 = l.get(6).toUpperCase();
                    seq2 = PStringUtils.fastSplitOnWhitespace(seq2).get(6).toUpperCase();
                    seq3 = PStringUtils.fastSplitOnWhitespace(seq3).get(6).toUpperCase();
                    for (int i = 0; i < seq1.length(); i++) {
                        sb.setLength(0);
                        char c = seq1.charAt(i);
                        if (c == '-') continue;
                        sb.append(chr).append("\t").append(pos).append("\t").append(c).append("\t");
                        sb.append(seq2.charAt(i)).append("\t").append(seq3.charAt(i));
                        pos++;
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            
        });
    }
}
