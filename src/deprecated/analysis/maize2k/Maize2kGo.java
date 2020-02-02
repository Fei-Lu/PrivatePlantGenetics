/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deprecated.analysis.maize2k;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.format.tree.Newick;


/**
 *
 * @author feilu
 */
public class Maize2kGo {
    
    public Maize2kGo () {
        //this.processReference();
        this.processAnnotation();
        //this.seqQualityTestPipe();
        //this.setHapScannerPipe();
        //this.test();
    }
    
    public void test () {
        String nwkS = "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);";
        Newick nwk = new Newick(nwkS);
    }
    
    private void setHapScannerPipe () {
        //new HapMap3Processor();
        //new HapMapTaxaProcessor();
        //String parameterFileS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/parameters_hapScanner.txt";
        String parameterFileS = "/Users/feilu/Documents/analysisL/production/parameters_hapScanner.txt";
        new HapScanner(parameterFileS);
    }
    
    private void seqQualityTestPipe () {
        new FastqQuality ();
    }
    
    private void processReference () {
        new ReferenceProcessor();
    }
    
    private void processAnnotation () {
        new AnnotationProcessor();
    }
    
    public static void main (String[] args) {
        new Maize2kGo();
    }
    
}
