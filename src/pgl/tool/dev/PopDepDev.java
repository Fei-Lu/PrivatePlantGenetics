package pgl.tool.dev;

import pgl.app.popdep.PopDep;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class PopDepDev {

    public PopDepDev () {
        //this.changeTaxaBamFormat();
//        this.runPopDep1();
        this.runPopDep2();
    }

    public void runPopDep2 () {
        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/parameters_popdep_2.txt";
        new PopDep(parameterFileS);
    }

    public void runPopDep1 () {
        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/parameters_popdep_1.txt";
        new PopDep(parameterFileS);
    }

    public void changeTaxaBamFormat () {
        String abInfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/source/004_taxaRefBam.ABgenome.removeBadTaxa.addNAFU.addS1.txt";
        String abdInfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/source/002_taxaRefBam.ABDgenome.manual.addNAFU.txt";
        String dInfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/source/003_taxaRefBam.Dgenome.addNAFU.txt";
        String abOutfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/source/vmap2_AB_taxaBam.txt";
        String abdOutfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/source/vmap2_ABD_taxaBam.txt";
        String dOutfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/source/vmap2_D_taxaBam.txt";
        this.changeFormat(abInfileS,abOutfileS);
        this.changeFormat(abdInfileS,abdOutfileS);
        this.changeFormat(dInfileS,dOutfileS);
    }

    private void changeFormat (String infileS, String outfileS) {
        RowTable<String> t = new RowTable<>(infileS);
        HashSet<String> set = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            set.add(t.getCell(i, 0));
        }
        String[] taxa = set.toArray(new String[set.size()]);
        Arrays.sort(taxa);
        List<String>[] bamList = new ArrayList[taxa.length];
        for (int i = 0; i < bamList.length; i++) {
            bamList[i] = new ArrayList<>();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(taxa, t.getCell(i, 0));
            bamList[index].add(t.getCell(i,2));
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tBams(A list of bams of the taxon, seperated by the delimiter of Tab)");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < taxa.length; i++) {
                sb.setLength(0);
                sb.append(taxa[i]);
                for (int j = 0; j < bamList[i].size(); j++) {
                    sb.append("\t").append(bamList[i].get(j));
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void commandDemo () {
        //samtools depth -r 1:30-40 -Q 20 /Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/bams/TW0063.sub.bam > /Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/a.txt
    }
}
