/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;
import utils.wheat.RefV1Utils;

/**
 *
 * @author feilu
 */
public class AncestralAnno {
    
    public AncestralAnno () {
        //this.splitMafBySubgenome();
        //this.splitMafByChr();
    }
    
    public void splitMafByChr () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/003_ancestral/bySubgenome/";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/003_ancestral/byChr";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".aln.gz");
        List<File> fList = Arrays.asList(fs);
        String header = "Chr\tPos\tAncestral\tRef\tUrartu\tBarley";
        String[] bases = {"A","G","T","C"};
        Arrays.sort(bases);
        fList.parallelStream().forEach(f -> {
            BufferedWriter[] bws = new BufferedWriter[14];
            String[] fileNames = new String[14];
            int[] chrIDs = new int[14];
            char subGenome = '-';
            if (f.getName().startsWith("A")) {
                subGenome = 'A';
                int cnt = 0;
                for (int i = 0; i < 7; i++) {
                    String fileName1 = "chr"+PStringUtils.getNDigitNumber(3, i*6+1)+"_ancestral.txt.gz";
                    fileName1 = new File(outDirS, fileName1).getAbsolutePath();
                    String fileName2 = "chr"+PStringUtils.getNDigitNumber(3, i*6+2)+"_ancestral.txt.gz";
                    fileName2 = new File(outDirS, fileName2).getAbsolutePath();
                    fileNames[cnt] = fileName1; 
                    chrIDs[cnt] = i*6+1;
                    cnt++;
                    fileNames[cnt] = fileName2; 
                    chrIDs[cnt] = i*6+2;
                    cnt++;
                }
            }
            else if (f.getName().startsWith("B")) {
                subGenome = 'B';
                int cnt = 0;
                for (int i = 0; i < 7; i++) {
                    String fileName1 = "chr"+PStringUtils.getNDigitNumber(3, i*6+3)+"_ancestral.txt.gz";
                    fileName1 = new File(outDirS, fileName1).getAbsolutePath();
                    String fileName2 = "chr"+PStringUtils.getNDigitNumber(3, i*6+4)+"_ancestral.txt.gz";
                    fileName2 = new File(outDirS, fileName2).getAbsolutePath();
                    fileNames[cnt] = fileName1; 
                    chrIDs[cnt] = i*6+3;
                    cnt++;
                    fileNames[cnt] = fileName2; 
                    chrIDs[cnt] = i*6+4;
                    cnt++;
                }
            }
            else if (f.getName().startsWith("D")) {
                subGenome = 'D';
                int cnt = 0;
                for (int i = 0; i < 7; i++) {
                    String fileName1 = "chr"+PStringUtils.getNDigitNumber(3, i*6+5)+"_ancestral.txt.gz";
                    fileName1 = new File(outDirS, fileName1).getAbsolutePath();
                    String fileName2 = "chr"+PStringUtils.getNDigitNumber(3, i*6+6)+"_ancestral.txt.gz";
                    fileName2 = new File(outDirS, fileName2).getAbsolutePath();
                    fileNames[cnt] = fileName1; 
                    chrIDs[cnt] = i*6+5;
                    cnt++;
                    fileNames[cnt] = fileName2; 
                    chrIDs[cnt] = i*6+6;
                    cnt++;
                }
            }
            try {
                for (int i = 0; i < bws.length; i++) {
                    bws[i] = IOUtils.getTextGzipWriter(fileNames[i]);
                    bws[i].write(header);
                    bws[i].newLine();
                }
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                StringBuilder sb = new StringBuilder();
                List<String> l = new ArrayList();
                int index = -1;
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(temp);
                    String chromosome = l.get(0).replaceFirst("chr", "")+subGenome;
                    int chrID = RefV1Utils.getChrID(chromosome, Integer.parseInt(l.get(1)));
                    int pos = RefV1Utils.getPosOnChrID(chromosome, Integer.parseInt(l.get(1)));
                    index = Arrays.binarySearch(chrIDs, chrID);
                    sb.append(chrID).append("\t").append(pos).append("\t");
                    if (l.get(3).equals(l.get(4)) && Arrays.binarySearch(bases, l.get(3)) > -1) {
                        sb.append(l.get(3)).append("\t");
                    }
                    else sb.append("NA").append("\t");
                    sb.append(l.get(2)).append("\t").append(l.get(3)).append("\t").append(l.get(4));
                    bws[index].write(sb.toString());
                    bws[index].newLine();
                }
                for (int i = 0; i < bws.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void splitMafBySubgenome () {
        String inDirS = "/Volumes/Fei_HDD_Mac/Gerp/asAlle/";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/003_ancestral/bySubgenome/";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".maf");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = f.getName().replaceFirst(".maf", ".aln.gz");
            outfileS = new File(outDirS, outfileS).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
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
