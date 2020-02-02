package analysis.wheatHapMap;

import pgl.format.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.graphcis.r.BoxPlot;
import pgl.utils.Dyad;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;
import pgl.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.LongAdder;

public class DeleteriousSite {

    public DeleteriousSite () {
//        this.mkSynSite();
//        this.mkNonSynSite();
//        this.mkNonSynSiftRefAltSite();
//        this.mkNonSynSiftDerivedSite();
        //this.mkGerpSite();
//        this.mkBurdenFiles();
//        this.plotBurden();
    }

    public void plotBurden () {
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite";
        List<File> fList = IOUtils.getFileListInDirEndsWith(dirS, ".txt");
        fList.parallelStream().forEach(f -> {
            String outpdf1 = f.getAbsolutePath().replaceFirst(".txt", "_count.pdf");
            String outpdf2 = f.getAbsolutePath().replaceFirst(".txt", "_countPerSite.pdf");
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            HashSet<String> typeSet = new HashSet<>();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append(t.getCell(i, 0)).append("_").append(t.getCell(i, 1));
                typeSet.add(sb.toString());
            }
            String[] types = typeSet.toArray(new String[typeSet.size()]);
            Arrays.sort(types);
            TDoubleArrayList[] alleleCountLists = new TDoubleArrayList[types.length];
            TDoubleArrayList[] ratioLists = new TDoubleArrayList[types.length];
            for (int i = 0; i < types.length; i++) {
                alleleCountLists[i] = new TDoubleArrayList();
                ratioLists[i] = new TDoubleArrayList();
            }
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);;
                sb.append(t.getCell(i, 0)).append("_").append(t.getCell(i, 1));
                int index = Arrays.binarySearch(types, sb.toString());
                alleleCountLists[index].add(Double.parseDouble(t.getCell(i, 3)));
                ratioLists[index].add(Double.parseDouble(t.getCell(i, 5)));
            }
            double[][] counts = new double[types.length][];
            double[][] ratios = new double[types.length][];
            for (int i = 0; i < types.length; i++) {
                counts[i] = alleleCountLists[i].toArray();
                ratios[i] = ratioLists[i].toArray();
            }
            BoxPlot plot = new BoxPlot(counts, types);
            plot.setTitle(f.getName()+"_alleleCount");
            plot.saveGraph(outpdf1);
            plot = new BoxPlot(ratios, types);
            plot.setTitle(f.getName()+"_alleleCountPerSite");
            plot.saveGraph(outpdf2);
            System.out.println(f.getName());
        });

    }

    public void mkBurdenFiles () {
        String vcfDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/002_exonSNPVCF";
        String taxaFileS = "/Users/feilu/Documents/analysisH/vmap2/006_populationStructure/vmap2_taxa.txt";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite";
        File[] fs = new File(dirS).listFiles();
        List<File> dirList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (!fs[i].isDirectory()) continue;
            dirList.add(fs[i]);
        }
        Collections.sort(dirList);
        RowTable<String> t = new RowTable<>(taxaFileS);
        int groupIndex = t.getColumnIndex("TreeValidatedGroupbyPloidy");
        List<String> gList = t.getColumn(groupIndex);
        HashSet<String> groupSet = new HashSet<>(gList);
        groupSet.remove("ExclusionHexaploid");
        groupSet.remove("ExclusionTetraploid");
        String[] groups = groupSet.toArray(new String[groupSet.size()]);
        Arrays.sort(groups);
        List<String>[] groupTaxaLists = new ArrayList[groups.length];
        String[] groupPloidy = new String[groups.length];
        for (int i = 0; i < groupTaxaLists.length; i++) {
            groupTaxaLists[i] = new ArrayList<>();
        }
        int index = -1;
        String query = null;
        for (int i = 0; i < t.getRowNumber(); i++) {
            query = t.getCell(i, groupIndex);
            index = Arrays.binarySearch(groups, query);
            if (index < 0) continue;
            groupTaxaLists[index].add(t.getCell(i, 4));
            groupPloidy[index] = t.getCell(i, 5);
        }
        for (int i = 0; i < groupTaxaLists.length; i++) {
            Collections.sort(groupTaxaLists[i]);
        }
        dirList.parallelStream().forEach(dir -> {
            this.outputBurden(dir, vcfDirS, groups, groupTaxaLists, groupPloidy, dirS);
        });

    }

    private void outputBurden (File delDir, String vcfDirS, String[] groups, List<String>[] groupTaxaLists, String[] groupPloidy, String outputDirS) {
        String outFileS = new File (outputDirS, delDir.getName()+".txt").getAbsolutePath();
        String[] genomeTypes = {"A", "B", "D"};

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outFileS);
            String header = "Subgenome\tPloidy\tTaxa\tDerivedAlleleCount\tScoredGenotypeCount\tDerivedAllelePerSite";
            bw.write(header);
            bw.newLine();
            int[] groupIndex = new int[2];
            for (int i = 0; i < genomeTypes.length; i++) {
                int cnt = 0;
                for (int j = 0; j < groupPloidy.length; j++) {
                    if (groupPloidy[j].contains(genomeTypes[i])) {
                        groupIndex[cnt] = j;
                        cnt++;
                    }
                }
                int[][] alleleCount = new int[groupIndex.length][];
                int[][] genotypeCount = new int[groupIndex.length][];
                int[][] taxaIndex = new int[groupIndex.length][];
                for (int j = 0; j < groupIndex.length; j++) {
                    alleleCount[j] = new int[groupTaxaLists[groupIndex[j]].size()];
                    genotypeCount[j] = new int[groupTaxaLists[groupIndex[j]].size()];
                    taxaIndex[j] = new int[groupTaxaLists[groupIndex[j]].size()];
                }
                int[] chrIDs = RefV1Utils.getChrIDsOfSubgenome(genomeTypes[i]);
                for (int j = 0; j < chrIDs.length; j++) {
                    String vcfFileS = "chr"+PStringUtils.getNDigitNumber(3, chrIDs[j])+"_exon_vmap2.1.vcf.gz";
                    vcfFileS = new File (vcfDirS, vcfFileS).getAbsolutePath();
                    String delFileS = "chr"+PStringUtils.getNDigitNumber(3, chrIDs[j])+"_SNP_anno.txt.gz";
                    delFileS = new File (delDir, delFileS).getAbsolutePath();
                    String temp = null;
                    List<String> l = new ArrayList<>();
                    TIntArrayList posList = new TIntArrayList();
                    List<String> ancestraList = new ArrayList<>();
                    BufferedReader br = IOUtils.getTextGzipReader(delFileS);
                    br.readLine();
                    String ances = null;
                    while((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp);
                        ances = l.get(15);
                        if (ances.startsWith("N")) continue;
                        posList.add(Integer.parseInt(l.get(2)));
                        ancestraList.add(ances);
                    }
                    br.close();
                    int[] delPos = posList.toArray();
                    String[] ancestral = ancestraList.toArray(new String[ancestraList.size()]);
                    br = IOUtils.getTextGzipReader(vcfFileS);
                    while ((temp = br.readLine()).startsWith("##")) {}
                    l = PStringUtils.fastSplit(temp);
                    int index = -1;
                    for (int k = 0; k < groupIndex.length; k++) {
                        for (int m = 0; m < l.size(); m++) {
                            index = Collections.binarySearch(groupTaxaLists[groupIndex[k]], l.get(m));
                            if (index < 0) continue;
                            taxaIndex[k][index] = m;
                        }
                    }
                    int derivedCode = -1;
                    String currentVCF = null;
                    String[] tem = null;
                    while ((temp = br.readLine()) != null) {
                        l = PStringUtils.fastSplit(temp.substring(0,60));
                        index = Arrays.binarySearch(delPos, Integer.parseInt(l.get(1)));
                        if (index < 0) continue;
                        l = PStringUtils.fastSplit(temp);
                        if (l.get(3).equals(ancestral[index])) derivedCode = 1;
                        else derivedCode = 0;
                        for (int k = 0; k < groupIndex.length; k++) {
                            for (int m = 0; m < groupTaxaLists[groupIndex[k]].size(); m++) {
                                currentVCF = l.get(taxaIndex[k][m]);
                                if (currentVCF.startsWith(".")) {
                                    continue;
                                }
                                genotypeCount[k][m]++;
                                tem = currentVCF.split(":")[0].split("/");
                                for (int n = 0; n < tem.length; n++) {
                                    if (Integer.parseInt(tem[n]) == derivedCode) {
                                        alleleCount[k][m]++;
                                    }
                                }
                            }
                        }
                    }
                    br.close();
                    System.out.println(vcfFileS);
                }
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < groupIndex.length; j++) {
                    for (int k = 0; k < alleleCount[j].length; k++) {
                        sb.setLength(0);
                        sb.append(genomeTypes[i]).append("\t").append(groupPloidy[groupIndex[j]]).append("\t");
                        sb.append(groupTaxaLists[groupIndex[j]].get(k)).append("\t");
                        sb.append(alleleCount[j][k]).append("\t").append(genotypeCount[j][k]).append("\t");
                        sb.append((float)((double)alleleCount[j][k]/genotypeCount[j][k]));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();

        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkGerpSite () {
        String exonSNPDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite/005_gerp";
        List<File> fList = IOUtils.getFileListInDirEndsWith(exonSNPDirS, "gz");
        double gerpThresh = 1;
        LongAdder counter = new LongAdder();
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String header = two.getFirstElement();
            List<String> records = two.getSecondElement();
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(header);
                bw.newLine();
                List<String> l = new ArrayList<>();
                int gerpIndex = VMapDBUtils.getColumnIndexInHeader(header, "Gerp");
                int cnt = 0;
                String derivedSift = null;
                for (int i = 0; i < records.size(); i++) {
                    l = PStringUtils.fastSplit(records.get(i));
                    if (l.get(gerpIndex).startsWith("NA")) continue;
                    bw.write(records.get(i));
                    bw.newLine();
                    cnt++;
                }
                bw.flush();
                bw.close();
                counter.add(cnt);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Total: " + String.valueOf(counter.longValue()));
    }

    public void mkNonSynSiftDerivedSite () {
        String exonSNPDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite/004_nonsynSiftDerived";
        List<File> fList = IOUtils.getFileListInDirEndsWith(exonSNPDirS, "gz");
        double siftThresh = 0.05;
        LongAdder counter = new LongAdder();
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String header = two.getFirstElement();
            List<String> records = two.getSecondElement();
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(header);
                bw.newLine();
                List<String> l = new ArrayList<>();
                int typeIndex = VMapDBUtils.getColumnIndexInHeader(header, "Variant_type");
                int derivedIndex = VMapDBUtils.getColumnIndexInHeader(header, "Derived_SIFT");
                int cnt = 0;
                String derivedSift = null;
                for (int i = 0; i < records.size(); i++) {
                    l = PStringUtils.fastSplit(records.get(i));
                    if (!l.get(typeIndex).startsWith("NONSYN")) continue;
                    derivedSift = l.get(derivedIndex);
                    if (derivedSift.startsWith("NA")) continue;
                    if (Double.parseDouble(derivedSift) < siftThresh) {
                        bw.write(records.get(i));
                        bw.newLine();
                        cnt++;
                    }
                }
                bw.flush();
                bw.close();
                counter.add(cnt);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Total: " + String.valueOf(counter.longValue()));
    }

    public void mkNonSynSiftRefAltSite () {
        String exonSNPDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite/003_nonsynSiftRefAlt";
        List<File> fList = IOUtils.getFileListInDirEndsWith(exonSNPDirS, "gz");
        double siftThresh = 0.05;
        LongAdder counter = new LongAdder();
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String header = two.getFirstElement();
            List<String> records = two.getSecondElement();
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(header);
                bw.newLine();
                List<String> l = new ArrayList<>();
                int typeIndex = VMapDBUtils.getColumnIndexInHeader(header, "Variant_type");
                int altIndex = VMapDBUtils.getColumnIndexInHeader(header, "Alt_SIFT");
                int refIndex = VMapDBUtils.getColumnIndexInHeader(header, "Ref_SIFT");
                int cnt = 0;
                String altSift = null;
                String refSift = null;
                for (int i = 0; i < records.size(); i++) {
                    l = PStringUtils.fastSplit(records.get(i));
                    if (!l.get(typeIndex).startsWith("NONSYN")) continue;
                    altSift = l.get(altIndex);
                    refSift = l.get(refIndex);
                    if (altSift.startsWith("NA") && refSift.startsWith("NA")) continue;
                    if (!altSift.startsWith("NA") && Double.parseDouble(altSift) < siftThresh) {
                        bw.write(records.get(i));
                        bw.newLine();
                        cnt++;
                    }
                    else if (!refSift.startsWith("NA") && Double.parseDouble(refSift) < siftThresh) {
                        bw.write(records.get(i));
                        bw.newLine();
                        cnt++;
                    }
                }
                bw.flush();
                bw.close();
                counter.add(cnt);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Total: " + String.valueOf(counter.longValue()));
    }

    public void mkNonSynSite () {
        String exonSNPDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite/002_nonsyn/";
        List<File> fList = IOUtils.getFileListInDirEndsWith(exonSNPDirS, "gz");
        LongAdder counter = new LongAdder();
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String header = two.getFirstElement();
            List<String> records = two.getSecondElement();
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(header);
                bw.newLine();
                List<String> l = new ArrayList<>();
                int typeIndex = VMapDBUtils.getColumnIndexInHeader(header, "Variant_type");
                int cnt = 0;
                for (int i = 0; i < records.size(); i++) {
                    l = PStringUtils.fastSplit(records.get(i));
                    if (!l.get(typeIndex).startsWith("NONSYN")) continue;
                    bw.write(records.get(i));
                    bw.newLine();
                    cnt++;
                }
                bw.flush();
                bw.close();
                counter.add(cnt);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Total: " + String.valueOf(counter.longValue()));
    }

    public void mkSynSite () {
        String exonSNPDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/004_deleteriousISite/001_syn";
        List<File> fList = IOUtils.getFileListInDirEndsWith(exonSNPDirS, "gz");
        LongAdder counter = new LongAdder();
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String header = two.getFirstElement();
            List<String> records = two.getSecondElement();
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(header);
                bw.newLine();
                List<String> l = new ArrayList<>();
                int typeIndex = VMapDBUtils.getColumnIndexInHeader(header, "Variant_type");
                int cnt = 0;
                for (int i = 0; i < records.size(); i++) {
                    l = PStringUtils.fastSplit(records.get(i));
                    if (!l.get(typeIndex).startsWith("SYN")) continue;
                    bw.write(records.get(i));
                    bw.newLine();
                    cnt++;
                }
                bw.flush();
                bw.close();
                counter.add(cnt);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Total: " + String.valueOf(counter.longValue()));
    }
}
