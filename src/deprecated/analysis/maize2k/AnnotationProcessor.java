/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package deprecated.analysis.maize2k;

import pgl.format.genomeAnnotation.GFFUtils;
import pgl.format.genomeAnnotation.GeneFeature;
import java.io.BufferedReader;
import java.util.HashSet;
import java.util.Set;
import pgl.utils.IOUtils;

/**
 *
 * @author feilu
 */
public class AnnotationProcessor {
    
    public AnnotationProcessor () {
        //this.modifyGFF3();
        //this.testField();
        this.mkPGFFile();
        //this.mkGeneAnnotationFileS();
    }
    
    private void mkGeneAnnotationFileS () {
        String inputGFF3 = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.modified.gff3.gz";
        String outputFileS = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.modified.geneAnnotation.txt";
        GFFUtils.mkMaizeGeneAnnotationFile(inputGFF3, outputFileS);
    }
    
    private void mkPGFFile () {
        String inputGFF3 = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.modified.gff3";
        String outputPGF = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.pgf";
        GeneFeature gf = new GeneFeature();
        gf.readFromMaizeGFF(inputGFF3);
        gf.writeFile(outputPGF);
    }
    
    private void testField () {
        String inputGFF3 = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.modified.gff3";
        Set<String> geneTypeSet = new HashSet<>();
        Set<String> bioTypeSet = new HashSet<>();
        String temp = null;
        try {
            BufferedReader br = IOUtils.getTextReader(inputGFF3);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("###")) {
                    temp = br.readLine();
                    if (temp == null) break;
                    String[] tem = temp.split("\t");
                    geneTypeSet.add(tem[2]);
                    tem = temp.split(";");
                    String[] te = temp.split("\t");
                    if (te[2].startsWith("gene") || te[2].startsWith("ncRNA_gene")) {
                        String biotype = "NA";
                        for (int i = 1; i < tem.length; i++) {
                            if (tem[i].startsWith("biotype")) biotype = tem[i].replaceFirst("biotype=", "");      
                            bioTypeSet.add(biotype);
                        }
                    }
                    
                }                   
            }
            br.close();
        }
        catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
        System.out.println("Genetype Set");
        for (String s : geneTypeSet) {
            System.out.println(s);
        }
        System.out.println("Biotype Set");
        for (String s : bioTypeSet) {
            System.out.println(s);
        }
    }
    
    private void modifyGFF3 () {
        String inputGFF3 = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.gff3.gz";
        String outputGFF3 = "/Users/feilu/Documents/database/maize/gene/Zea_mays.AGPv4.38.modified.gff3";
        GFFUtils.modifyMaizeAGPV4GFF3(inputGFF3, outputGFF3);
    }
}
