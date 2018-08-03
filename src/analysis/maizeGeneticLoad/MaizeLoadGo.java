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
        this.mkHmp3Header();
        this.processGERP();
    }
    
    public void mkHmp3Header () {
        
    }
    
    public void processGERP () {
        new GERP();
    }
    
    public static void main (String[] args) {
        new MaizeLoadGo ();
    }
    
}
