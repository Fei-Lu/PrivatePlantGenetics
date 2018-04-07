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
        this.test();
    }
    
    public void test () {
        
        }
    
    public void referenceGenome () {
        new WheatReferenceGenome();
    }
    
    public static void main (String[] args) {
        new WheatHmpGo();
    }
}
