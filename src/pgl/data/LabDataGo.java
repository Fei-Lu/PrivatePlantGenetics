/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.data;

/**
 *
 * @author feilu
 */
public class LabDataGo {
    
    public LabDataGo() {
        //this.referenceGenome();
//        this.referenceGene();
//        this.bamManagement();
    }

    public void bamManagement () {
        new BAMManagement();
    }

    public void referenceGene () {
        new WheatReferenceGene();
    }
    
    public void referenceGenome () {
        new WheatReferenceGenome();
    }


}
