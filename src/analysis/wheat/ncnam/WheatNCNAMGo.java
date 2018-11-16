/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.ncnam;

/**
 *
 * @author feilu
 */
class WheatNCNAMGo {
    
    public WheatNCNAMGo () {
        this.selectTaxa();
    }
    
    public void selectTaxa () {
        new TaxaSelection();
    }
    
    public static void main (String[] args) {
        new WheatNCNAMGo ();
    }
    
}
