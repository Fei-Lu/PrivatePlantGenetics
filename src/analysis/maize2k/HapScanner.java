/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maize2k;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import zhouyao.analysis.wheatHapMap.IOUtils;

/**
 *
 * @author feilu
 */
public class HapScanner {
    //The path of reference genome with .fai file in the same folder
    String referenceGenomeS = null;
    //The directory of sorted bam files, with .bai files in the same folder
    String bamDirS = null;
    //The directory of posAllele files (with header), the format is Chr\tPos\tRef\tAlt (from VCF format). The positions come from haplotype library.
    String posAlleleDirS = null;
    //The directory of pos files (without header), the format is Chr\tPos. The positions come from haplotype library.
    String posDirS = null;
    //The chromosome which will be genotyped
    int chr = -1;
    //The path of samtools
    String samtoolsPath = null;
    //The directory of output
    String outputDirS = null;
    
    public HapScanner (String infileS) {
        this.parseParameters(infileS);
    }
    
    public void parseParameters (String infileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("HapScanner")) ifOut = true;
            if (!(temp = br.readLine()).equals("Author: Aoyue Bi, Fei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Email: biaoyue17@genetics.ac.cn; flu@genetics.ac.cn")) ifOut = true;
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) ifOut = true;
            if (ifOut) {
                System.out.println("Thanks for using HapScanner.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        referenceGenomeS = pLineList.get(0);
        bamDirS = pLineList.get(1);
        posAlleleDirS = pLineList.get(2);
        posDirS = pLineList.get(3);
        chr = Integer.valueOf(pLineList.get(4));
        samtoolsPath = pLineList.get(5);
        outputDirS = pLineList.get(6);   
    }
    
}
