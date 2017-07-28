/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maize2k;

import format.Fasta;
import java.io.File;
import java.util.ArrayList;

/**
 *
 * @author feilu
 */
public class ReferenceProcessor {
    
    public ReferenceProcessor () {
        //this.merge();
        //this.statistics();
        
    }
    
    private void merge () {
        String inputDirS = "/Users/feilu/Documents/database/maize/reference/download";
        String outputFileS = "/Users/feilu/Documents/database/maize/reference/merged/maizeAGPv4.fa";
        ArrayList<Fasta> faList = new ArrayList();
        String chr0FileS =  "/Users/feilu/Documents/database/maize/reference/download/Zea_mays.AGPv4.dna.nonchromosomal.fa.gz";
        Fasta fa = new Fasta (chr0FileS);
        faList.add(fa);
        for (int i = 0; i < 10; i++) {
            String inputFileS = new File (inputDirS, "Zea_mays.AGPv4.dna.chromosome." +String.valueOf(i+1)+ ".fa.gz").getAbsolutePath();
            fa = new Fasta (inputFileS);
            faList.add(fa);
        }
        String mtFileS =  "/Users/feilu/Documents/database/maize/reference/download/Zea_mays.AGPv4.dna.chromosome.Mt.fa.gz";
        fa = new Fasta (mtFileS);
        faList.add(fa);
        String ptFileS =  "/Users/feilu/Documents/database/maize/reference/download/Zea_mays.AGPv4.dna.chromosome.Pt.fa.gz";
        fa = new Fasta (ptFileS);
        faList.add(fa);
        String[] names = new String[faList.size()];
        String[] seqs = new String[faList.size()];
        int[] ids = new int[faList.size()];
        String ns = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        Fasta[] fas = faList.toArray(new Fasta[faList.size()]);
        for (int i = 0; i < fas.length; i++) {
            names[i] = String.valueOf(i);
            ids[i] = i + 1;
            StringBuilder sb = new StringBuilder(fas[i].getSeq(0));
            for (int j = 0; j < fas[i].getSeqNumber()-1; j++) {
                sb.append(ns).append(fas[i].getSeq(j+1));
            }
            seqs[i] = sb.toString();
        }
        fa = new Fasta (names, seqs, ids);
        System.out.println(fa.getSeqLength(0));
        fa.writeFasta(outputFileS);
    }
    
    private void statistics () {
        String infileS = "/Users/feilu/Documents/database/maize/reference/merged/maizeAGPv4.fa";
        Fasta f = new Fasta(infileS);
        System.out.println("Chromosomes\tLength");
        for (int i = 0; i < f.getSeqNumber(); i++) {
            System.out.println(String.valueOf(i)+"\t"+String.valueOf(f.getSeqLength(i)));
        }
    }
    
}
