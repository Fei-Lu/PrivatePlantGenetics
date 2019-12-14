/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import java.io.File;

/**
 *
 * @author feilu
 */
public class DeleteriousBiology {
    
    public DeleteriousBiology () {
        this.mkVCF();
    }
    
    public void mkVCF () {
        double siftThresh = 0.05;
        double gerpThresh = 1;
        double phyloPThresh = 0.5;
        String indirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/genicSNPAnnotation";
        String[] fs = new File(indirS).list();
        
    }
}
