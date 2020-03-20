/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.VMap2;

/**
 *
 * @author feilu
 */
public class WheatVMap2Go {
    
    public WheatVMap2Go() {
//        this.VMapII();
//        this.DBGene();
//        this.DBDeleterious();
//        this.DBWindow();
//        this.annotation();
//        this.deleteriousSite();
//        this.populationStructure();
        //this.deleteriousBiology();

    }
    
    public void deleteriousBiology () {
        new DeleteriousBiology ();
    }

    public void deleteriousSite () {
        new DeleteriousSite ();
    }

    public void populationStructure () {
        new PopulationStructure();
    }

    public void annotation () {
        this.annoAncestral();
//        this.annoSift();
//        this.annoGerp();
//        this.annoPhyloP();
//        this.annoCrossover();

    }

    public void annoCrossover () {
        new AnnoCrossover();
    }

    public void annoPhyloP() {
        new AnnoPhyloP();
    }
    
    public void annoGerp () {
        new AnnoGerp();
    }
    
    public void annoSift() {
        new AnnoSift();
    }
    
    public void annoAncestral() {
        new AnnoAncestral();
    }

    public void DBWindow () {
        new DBWindow ();
    }

    public void DBDeleterious () {
        new DBDeleterious();
    }
    
    public void DBGene() {
        new DBGene();
    }
    
    public void VMapII () {
        new VMapII();
    }

}
