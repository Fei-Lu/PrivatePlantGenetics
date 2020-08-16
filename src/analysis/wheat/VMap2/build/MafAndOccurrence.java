package analysis.wheat.VMap2.build;

import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

class MafAndOccurrence {

    public MafAndOccurrence () {
        this.filter();
    }
    //didn't complete here, this was done by aoyue
    public void filter () {
        String inputVCFDirS = "/Volumes/VMap2_Fei/vcf/003_hapScanner/";
        String outputVCFDirS = "/Volumes/VMap2_Fei/vcf/003_hapScanner/output";
        String ABTaxaFileS = "/Volumes/VMap2_Fei/vcf/003_hapScanner/taxa/EmmerWheat_S187.txt";
        String ABDTaxaFileS = "/Volumes/VMap2_Fei/vcf/003_hapScanner/taxa/BreadWheat_S420.txt";
        String DTaxaFileS = "/Volumes/VMap2_Fei/vcf/003_hapScanner/taxa/Ae.tauschii_S36.txt";

        int occu = 2;
        float mafThresh = (float)0.01;
        float missingThresh = (float)0.2;
        List<File> fList = IOUtils.getFileListInDirEndsWith(inputVCFDirS, ".vcf");
        Collections.sort(fList);
        fList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append(f.getName().split("\\.")[0]).append("_occu").append(occu).append("_maf").append(mafThresh).append("_miss").append(missingThresh).append(".vcf");
            String outfileS = new File (outputVCFDirS, sb.toString()).getAbsolutePath();
            int chr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chr);
            String[] abTaxa = null;
            String[] abdTaxa = null;
            String[] dTaxa = null;
            List<String> tList = new ArrayList<>();
            if (subgenome.equals("A") || subgenome.equals("B")) {
                try {
                    tList = new ArrayList<>();
                    BufferedReader br = IOUtils.getTextReader(ABTaxaFileS);
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        tList.add(temp);
                    }
                    br.close();
                    abTaxa = tList.toArray(new String[tList.size()]);
                    Arrays.sort(abTaxa);
                    tList = new ArrayList<>();
                    br = IOUtils.getTextReader(ABDTaxaFileS);
                    while ((temp = br.readLine()) != null) {
                        tList.add(temp);
                    }
                    br.close();
                    abdTaxa = tList.toArray(new String[tList.size()]);
                    Arrays.sort(abdTaxa);
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            else if (subgenome.equals("D")){
                try {
                    tList = new ArrayList<>();
                    BufferedReader br = IOUtils.getTextReader(DTaxaFileS);
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        tList.add(temp);
                    }
                    br.close();
                    dTaxa = tList.toArray(new String[tList.size()]);
                    Arrays.sort(dTaxa);
                    tList = new ArrayList<>();
                    br = IOUtils.getTextReader(ABDTaxaFileS);
                    while ((temp = br.readLine()) != null) {
                        tList.add(temp);
                    }
                    br.close();
                    abdTaxa = tList.toArray(new String[tList.size()]);
                    Arrays.sort(abdTaxa);
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            GenotypeGrid gt = new GenotypeGrid(f.getAbsolutePath(), GenoIOFormat.VCF);
            gt.sortByTaxa();
            String[] allTaxa = gt.getTaxaNames();
            int[][] taxaIndices = new int[2][];
            String[][] subTaxa = new String[2][];
            GenotypeGrid[] gts = new GenotypeGrid[2];
            if (subgenome.equals("D")) {
                taxaIndices[0] = new int[abdTaxa.length];
                taxaIndices[1] = new int[dTaxa.length];
                subTaxa[0] = abdTaxa;
                subTaxa[1] = dTaxa;
            }
            else {
                taxaIndices[0] = new int[abTaxa.length];
                taxaIndices[1] = new int[abdTaxa.length];
                subTaxa[0] = abTaxa;
                subTaxa[1] = abdTaxa;
            }
            for (int i = 0; i < taxaIndices.length; i++) {
                for (int j = 0; j < taxaIndices[i].length; j++) {
                    taxaIndices[i][j] = Arrays.binarySearch(allTaxa, subTaxa[i][j]);
                    if (taxaIndices[i][j] < 0) System.out.println(subTaxa[i][j]);
                }
            }
            for (int i = 0; i < gts.length; i++) {
                gts[i] = GenotypeOperation.getSubsetGenotypeByTaxon(gt, taxaIndices[i]);
            }
            TIntArrayList posList = new TIntArrayList();
            for (int i = 0; i < gt.getSiteNumber(); i++) {
                double[] missing = new double[2];
                for (int j = 0; j < 2; j++) {
                    missing[j] = (double)gts[j].getMissingNumberBySite(i)/gts[j].getTaxaNumber();
                }
                if (missing[0] > missingThresh && missing[1] > missingThresh) continue;
                for (int j = 0; j < gts.length; j++) {

                    double maf = gts[j].getMinorAlleleFrequency(i);
                    if (!(maf < mafThresh)) {
                        posList.add(gt.getPosition(i));
                        break;
                    }
                    else if (!(gts[j].getAlternativeAlleleNumberBySite(i) < occu)) {
                        posList.add(gt.getPosition(i));
                        break;
                    }
                }
            }
            int[] positions = posList.toArray();
            System.out.println(positions.length);




        });

    }

}
