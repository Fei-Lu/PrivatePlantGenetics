/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maize2k;

/**
 *
 * @author feilu
 */
public class Maize2kGo {
    
    public Maize2kGo () {
        this.seqQualityTestPipe();
        //this.processReference();
    }
    
    
    private void seqQualityTestPipe () {
        new FastqQuality ();
    }
    
    private void processReference () {
        new ReferenceProcessor();
    }
    
    public static void main (String[] args) {
        new Maize2kGo();
    }
    
}
