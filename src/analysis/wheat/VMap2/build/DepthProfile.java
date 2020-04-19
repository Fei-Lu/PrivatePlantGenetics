package analysis.wheat.VMap2.build;

import it.unimi.dsi.fastutil.Hash;
import pgl.app.popdep.PopDep;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import sun.util.resources.cldr.zh.CalendarData_zh_Hans_HK;

import java.io.BufferedWriter;
import java.util.*;

class DepthProfile {

    public DepthProfile () {
        //this.mkTaxaBamFile();
        //this.mkChrLengthFile();
        this.runStep1();
    }

    public void runStep1 () {
        String abdInfileS = "/data1/home/feilu/popdep/abd_parameters_popdep_1.txt";
        String abInfileS = "/data1/home/feilu/popdep/ab_parameters_popdep_1.txt";
        String dInfileS = "/data1/home/feilu/popdep/d_parameters_popdep_1.txt";
        new PopDep(dInfileS);
        new PopDep(abInfileS);
        new PopDep(abdInfileS);
    }

    public void mkChrLengthFile () {
        String abdLengthFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/abdChromosomeLength.txt";

        try {
            BufferedWriter bw = IOUtils.getTextWriter(abdLengthFileS);
            bw.write("Chromosome(Should be all choromosomes)\tLength");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            int[] ids = RefV1Utils.getChrIDs();
            for (int i = 0; i < RefV1Utils.getChrIDs().length; i++) {
                sb.setLength(0);
                sb.append(ids[i]).append("\t").append(RefV1Utils.getChrIDLength(ids[i]));
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

    public void mkTaxaBamFile () {
        String vmap2TaxaFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/wheatVMapII_germplasmInfo_20191225.txt";
        String old2newFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/bam_Old2New.map.txt";
        String abdFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/002_taxaRefBam.ABDgenome.manual.addNAFU.txt";
        String abFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/003_taxaRefBam.Dgenome.addNAFU.txt";
        String dFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/source/004_taxaRefBam.ABgenome.removeBadTaxa.addNAFU.addS1.txt";
        String abdOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/abdTaxaBam.txt";
        String abOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/abTaxaBam.txt";
        String dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step1/dTaxaBam.txt";
        RowTable<String> t = new RowTable<>(old2newFileS);
        HashMap<String, String> old2new = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCell(i,0).endsWith(".bai")) continue;
            old2new.put(t.getCell(i,0), t.getCell(i,1));
        }
        t = new RowTable<>(vmap2TaxaFileS);
        HashSet<String> set = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.getCell(i,25).equals("1")) continue;
            set.add(t.getCell(i, 4));
        }
        String[] taxa = set.toArray(new String[set.size()]);
        Arrays.sort(taxa);
        System.out.println(taxa.length);
        this.mkTaxaBam(old2new, taxa, abdFileS, abdOutfileS);
        this.mkTaxaBam(old2new, taxa, abFileS, abOutfileS);
        this.mkTaxaBam(old2new, taxa, dFileS, dOutfileS);

    }

    private void mkTaxaBam (HashMap<String, String> old2new, String[] taxa, String infileS, String outfileS) {
        RowTable<String> t = new RowTable<>(infileS);
        List<String> taxaList = t.getColumn(0);
        HashSet<String> set = new HashSet<>(taxaList);
        set.retainAll(Arrays.asList(taxa));
        String[] nTaxa = set.toArray(new String[set.size()]);
        Arrays.sort(nTaxa);
        List<String>[] bamLists = new ArrayList[nTaxa.length];
        for (int i = 0; i < nTaxa.length; i++) {
            bamLists[i] = new ArrayList<>();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(nTaxa, t.getCell(i, 0));
            if (index < 0) continue;
            bamLists[index].add(old2new.get(t.getCell(i,2)));
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa	Bams(A list of bams of the taxon, seperated by the delimiter of Tab)");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < nTaxa.length; i++) {
                sb.setLength(0);
                sb.append(nTaxa[i]);
                for (int j = 0; j < bamLists[i].size(); j++) {
                    sb.append("\t").append(bamLists[i].get(j));
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
}
