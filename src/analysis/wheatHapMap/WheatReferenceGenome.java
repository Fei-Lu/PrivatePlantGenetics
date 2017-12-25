/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import com.koloboke.collect.map.hash.HashByteByteMap;

import format.dna.DNAUtils;
import format.dna.Fasta;

import java.util.Arrays;
import utils.Benchmark;

/**
 *
 * @author feilu
 */
class WheatReferenceGenome {
    
    public WheatReferenceGenome () {
        this.checkIUPAC();
        //this.splitGenome();
        
    }
    
    public void splitGenome () {
        String infileS = "/Users/feilu/Documents/database/wheat/reference/download/unzipped/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta";
        
        long s2 = System.nanoTime();
        for (int j = 0; j < 10000000; j++) {
            byte[] a = infileS.toUpperCase().getBytes();
        }
        System.out.println("t:"+(System.nanoTime()-s2));
        
        HashByteByteMap map = DNAUtils.getBaseLowerToUpperMap();
        long s1 = System.nanoTime();
        for (int j = 0; j < 10000000; j++) {
            byte[] b = infileS.getBytes();
            for (int i = 0; i < b.length; i++) {
                if (b[i] < 97) continue;
                if (b[i] > 122) continue;
                b[i] = map.get(b[i]);
            }
        }
        System.out.println("m:"+(System.nanoTime()-s1));
        
        
    }
    
    public void checkIUPAC () {
        String infileS = "/Users/feilu/Documents/database/wheat/reference/download/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta";
        Fasta f = new Fasta(infileS);
        System.out.println(f.isThereNonACGTNBase());
        
    }
    
    
}
