/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maizeGeneticLoad;

import utils.CrossMapUtils;

/**
 *
 * @author feilu
 */
public class MaizeLoadGo {
    
    public MaizeLoadGo () {
        this.analysisPipe();
    }
    
    public void analysisPipe () {
        //this.annotateVariants();
        this.processGERP();
        //this.processSIFT();
    }
    
    public void annotateVariants () {
        new VariantsAnnotation ();
    }
    
    public void processGERP () {
        new GERP();
    }
    
    public void processSIFT () {
        new SIFT ();
    }
    
    public static void main (String[] args) {
        new MaizeLoadGo ();
    }
    
}
