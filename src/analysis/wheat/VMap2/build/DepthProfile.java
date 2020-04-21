package analysis.wheat.VMap2.build;

import gnu.trove.list.array.TIntArrayList;
import it.unimi.dsi.fastutil.Hash;
import pgl.app.popdep.PopDep;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import sun.util.resources.cldr.zh.CalendarData_zh_Hans_HK;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

class DepthProfile {

    public DepthProfile () {
        //this.mkTaxaBamFile();
        //this.mkChrLengthFile();
//        this.runStep1();
//        this.mkStep2FileS();
//        this.runStep2();
    }

    public void runStep2 () {
        String parameterDirS = "/data1/home/feilu/popdepStep2/parameters";
        List<File> fList = IOUtils.getFileListInDirEndsWith(parameterDirS, ".txt");
        for (int i = 0; i < fList.size(); i++) {
            new PopDep(fList.get(i).getAbsolutePath());
            System.out.println(String.valueOf(i+1)+" "+fList.get(i).getName());
        }
        System.out.println("Done");
    }

    public void mkStep2FileS() {
        String sourceFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popdep/inputfile/parameters_popdep_2.txt";
        String dTaxaModeFileS = "/data1/home/feilu/popdep/dTaxaDepthMode.txt";
        String abTaxaModeFileS = "/data1/home/feilu/popdep/abTaxaDepthMode.txt";
        String abdTaxaModeFileS = "/data1/home/feilu/popdep/abdTaxaDepthMode.txt";
        String dTaxaBamFileS = "/data1/home/feilu/popdep/dTaxaBam.txt";
        String abTaxaBamFileS = "/data1/home/feilu/popdep/abTaxaBam.txt";
        String abdTaxaBamFileS = "/data1/home/feilu/popdep/abdTaxaBam.txt";
        String samPath = "/data1/programs/samtools-1.8/samtools";
        String dOutDirS = "/data1/home/feilu/popdepStep2/D";
        String abOutDirS = "/data1/home/feilu/popdepStep2/AB";
        String abdOutDirS = "/data1/home/feilu/popdepStep2/ABD";
        String step2DirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/step2";
        List<String> template = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(sourceFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                template.add(temp);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.mkStep2(template, dTaxaBamFileS, dTaxaModeFileS, samPath, dOutDirS, step2DirS, "D");
        this.mkStep2(template, abTaxaBamFileS, abTaxaModeFileS, samPath, abOutDirS, step2DirS, "AB");
        this.mkStep2(template, abdTaxaBamFileS, abdTaxaModeFileS, samPath, abdOutDirS, step2DirS, "ABD");
    }

    private void mkStep2 (List<String> template, String taxaBamFileS, String taxaModeFileS, String samPath, String outDirS, String step2, String genomeType) {
        StringBuilder sb = new StringBuilder();
        TIntArrayList indexList = new TIntArrayList();
        for (int i = 0; i < template.size(); i++) {
            String current = template.get(i);
            if (current.startsWith("@")) continue;
            else if (current.startsWith("#")) continue;
            else if (current.equals("")) continue;
            indexList.add(i);
        }
        int[] indices = indexList.toArray(new int[indexList.size()]);
        int[] chrIDs = null;
        TIntArrayList chrIDList = new TIntArrayList();
        if (genomeType.equals("AB")) {
            chrIDList.addAll(RefV1Utils.getChrIDsOfSubgenomeA());
            chrIDList.addAll(RefV1Utils.getChrIDsOfSubgenomeB());
        }
        else if (genomeType.equals("ABD")) {
            chrIDList.addAll(RefV1Utils.getChrIDs());
        }
        else if (genomeType.equals("D")) {
            chrIDList.addAll(RefV1Utils.getChrIDsOfSubgenomeD());
        }
        chrIDs = chrIDList.toArray();
        Arrays.sort(chrIDs);
        for (int i = 0; i < chrIDs.length; i++) {
            sb.setLength(0);
            sb.append("chr").append(PStringUtils.getNDigitNumber(3, chrIDs[i])).append("_").append(genomeType).append("_popdep_vmap2.txt.gz");
            String outfileS = new File(outDirS,sb.toString()).getAbsolutePath();
            sb.setLength(0);
            sb.append(genomeType).append("_").append(chrIDs[i]).append("_parameters_popdep_2.txt");
            String step2FileS = new File(step2, sb.toString()).getAbsolutePath();
            try {
                int cnt = 0;
                BufferedWriter bw = IOUtils.getTextWriter(step2FileS);
                for (int j = 0; j < template.size(); j++) {
                    int index = Arrays.binarySearch(indices, j);
                    if (index < 0) {
                        bw.write(template.get(j));
                    }
                    else {
                        if (cnt == 0) {
                            bw.write(taxaBamFileS);
                        }
                        else if (cnt == 1) {
                            bw.write(taxaModeFileS);
                        }
                        else if (cnt == 2) {
                            bw.write(String.valueOf(chrIDs[i]));
                        }
                        else if (cnt == 3) {
                            bw.write(String.valueOf(RefV1Utils.getChrIDLength(chrIDs[i])));
                        }
                        else if (cnt == 4) {
                            bw.write("5");
                        }
                        else if (cnt == 5) {
                            bw.write(samPath);

                        }
                        else if (cnt == 6) {
                            bw.write("32");
                        }
                        else if (cnt == 7) {
                            bw.write(outfileS);
                        }
                        cnt++;
                    }
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

    public void runStep1 () {
        String abdInfileS = "/data1/home/feilu/popdep/abd_parameters_popdep_1.txt";
        String abInfileS = "/data1/home/feilu/popdep/ab_parameters_popdep_1.txt";
        String dInfileS = "/data1/home/feilu/popdep/d_parameters_popdep_1.txt";
        //new PopDep(dInfileS);
//        new PopDep(abInfileS);
//        new PopDep(abdInfileS);
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
