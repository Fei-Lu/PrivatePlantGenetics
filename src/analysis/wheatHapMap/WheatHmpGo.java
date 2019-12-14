/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import com.koloboke.collect.map.hash.HashByteByteMap;
import format.dna.BaseEncoder;
import format.table.ColumnTable;
import format.table.RowTable;
import format.table.TableInterface;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author feilu
 */
class WheatHmpGo {
    
    public WheatHmpGo () {
        //this.referenceGenome();
        //this.VMapII();
        //this.geneDB();
        //this.deleteriousDB();
        //this.annotation();
        //this.deleteriousBiology();
    }
    
    public void deleteriousBiology () {
        //new DeleteriousBiology ();
    }
    
    public void annotation () {
        //this.ancestralAnno();
        this.siftAnno();
        //this.gerpAnno();
        //this.phyloPAnno();
    }
    
    public void phyloPAnno () {
        new PhyloPAnno();
    }
    
    public void gerpAnno () {
        new GerpAnno ();
    }
    
    public void siftAnno () {
        new SiftAnno ();
    }
    
    public void ancestralAnno () {
        new AncestralAnno();
    }
    
    public void deleteriousDB () {
        new DeleteriousDB();
    }
    
    public void geneDB () {
        new GeneDB ();
    }
    
    public void VMapII () {
        new VMapII();
    }
    
    public static void main (String[] args) {
        new WheatHmpGo();
    }
}
