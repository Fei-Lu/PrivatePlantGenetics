/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.publicData;

/**
 *
 * @author feilu
 */
public class LabPublic {
    
    public LabPublic () {
        //this.referenceGenome();
        this.referenceGene();
    }
    
    public void referenceGene () {
        new WheatReferenceGene();
    }
    
    public void referenceGenome () {
        new WheatReferenceGenome();
    }
    
    public static void main (String[] args) {
        new LabPublic ();
    }
}
