/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.RNASeq;

/**
 *
 * @author feilu
 */
public class WheatRNASeqGo {
    
    public WheatRNASeqGo () {
        this.doMiscellaneous();
    }
    
    public void doMiscellaneous () {
        new Miscellaneous ();
    }
    
    public static void main (String[] args) {
        new WheatRNASeqGo ();
    }
}
