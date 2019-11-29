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
        this.ancestralDB();
    }
    
    public void ancestralDB () {
        new AncestralDB();
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
