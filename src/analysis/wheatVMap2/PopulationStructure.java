package analysis.wheatVMap2;

import pgl.format.table.RowTable;
import pgl.utils.IOFileFormat;
import pgl.utils.IOUtils;

import java.io.BufferedReader;
import java.util.Arrays;
import java.util.HashSet;

public class PopulationStructure {

    public PopulationStructure () {
        this.mkVMap2TaxaGroup();
    }

    public void mkVMap2TaxaGroup () {
        String infileS = "/Users/feilu/Documents/analysisH/vmap2/006_populationStructure/source/wheatVMapII_germplasmInfo_20191225.txt";
        String vcfFileS1 = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/002_exonSNPVCF/chr001_exon_vmap2.1.vcf.gz";
        String vcfFileS2 = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/002_exonSNPVCF/chr005_exon_vmap2.1.vcf.gz";
        String outfileS = "/Users/feilu/Documents/analysisH/vmap2/006_populationStructure/vmap2_taxa.txt";
        try {
            BufferedReader br = IOUtils.getTextGzipReader(vcfFileS1);
            String temp = null;
            while ((temp=br.readLine()).startsWith("##")) {}
            br.close();
            String[] tem = temp.split("\t");
            HashSet<String> vcfTaxaSet = new HashSet<>();
            for (int i = 0; i < tem.length-9; i++) {
                vcfTaxaSet.add(tem[i+9]);
            }
            br = IOUtils.getTextGzipReader(vcfFileS2);
            while ((temp=br.readLine()).startsWith("##")) {}
            br.close();
            tem = temp.split("\t");
            for (int i = 0; i < tem.length-9; i++) {
                vcfTaxaSet.add(tem[i+9]);
            }
            String[] vcfTaxa = vcfTaxaSet.toArray(new String[vcfTaxaSet.size()]);
            Arrays.sort(vcfTaxa);
            RowTable<String> t = new RowTable<>(infileS);
            int taxaIndex = t.getColumnIndex("Taxa");
            HashSet<String> taxaSet = new HashSet<>(t.getColumn(taxaIndex));
            boolean[] ifMarked = new boolean[vcfTaxa.length];
            boolean[] ifOut = new boolean[t.getRowNumber()];
            int index = -1;
            for (int i = 0; i < t.getRowNumber(); i++) {
                index = Arrays.binarySearch(vcfTaxa, t.getCell(i, taxaIndex));
                if (index < 0) continue;
                if (ifMarked[index]) continue;
                ifMarked[index] = true;
                ifOut[i] = true;
            }
            t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
            System.out.println("Total taxa: "+ String.valueOf(vcfTaxa.length));
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }
}
