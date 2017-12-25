/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

/**
 *
 * @author feilu
 */
class WheatHmpGo {
    
    public WheatHmpGo () {
        this.referenceGenome();
    }
    
    public void referenceGenome () {
        new WheatReferenceGenome();
    }
    
    public static void main (String[] args) {
        new WheatHmpGo();
    }
}
