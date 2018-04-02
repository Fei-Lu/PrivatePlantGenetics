/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maize2k;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IOUtils;


/**
 *
 * @author feilu
 */
public class Maize2kGo {
    
    public Maize2kGo () {
        //this.processReference();
        //this.processAnnotation();
        //this.seqQualityTestPipe();
        this.setHapScannerPipe();
        //this.test();
    }
    
    private void setHapScannerPipe () {
        //new HapMap3Processor();
        new HapMapTaxaProcessor();
//        String parameterFileS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/parameters_hapScanner.txt";
//        new HapScanner(parameterFileS);
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
