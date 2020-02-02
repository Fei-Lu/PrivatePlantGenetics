/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatVMap2;

import pgl.format.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import pgl.utils.IOUtils;
import pgl.utils.PArrayUtils;
import pgl.utils.PStringUtils;

/**
 *
 * @author feilu
 */
public class VMapII {

    public VMapII() {
        //this.mergePosList();
        this.mergeVCFandFilter();
    }
    
    //not done yet
    public void mergeVCFandFilter() {
        double mafThresh = 0.01;
        double missingThresh = 0.2;
        int vcfFileNumber = 2;
        String taxaFile1 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergeVCFfilter/EmmerWheat_S187.txt";
        String taxaFile2 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergeVCFfilter/BreadWheat_S419.txt";
        String vcfFile1 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergeVCFfilter/chr001_tetraploid.vcf";
        String vcfFile2 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergeVCFfilter/chr001_hexaploid.vcf";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergeVCFfilter/chr001_vmap2.vcf";
        List<String>[] taxaLists = new ArrayList[vcfFileNumber];
        List<String>[] taxaListOlds = new ArrayList[vcfFileNumber];
        RowTable[] ts = new RowTable[vcfFileNumber];
        ts[0] = new RowTable(taxaFile1);
        ts[1] = new RowTable(taxaFile2);
        int[][] columnIndices = new int[vcfFileNumber][];
        List<String> allTaxaList = new ArrayList();
        for (int i = 0; i < vcfFileNumber; i++) {
            ts[i].sortAsText(0);
            taxaLists[i] = ts[i].getColumn(0);
            taxaListOlds[i] = ts[i].getColumn(1);
            allTaxaList.addAll(taxaLists[i]);
            columnIndices[i] = new int[taxaLists[i].size()];
        }
//        ts[0].writeTextTable(taxaFile1, IOFileFormat.Text);
//        ts[1].writeTextTable(taxaFile2, IOFileFormat.Text);
        Collections.sort(allTaxaList);
        int[] tableIndex = new int[allTaxaList.size()];
        int[] columnIndex = new int[allTaxaList.size()];
        for (int i = 0; i < tableIndex.length; i++) {
            int index = Collections.binarySearch(taxaLists[0], allTaxaList.get(i));
            if (index < 0) {
                tableIndex[i] = 1;
            }
        }
        BufferedReader[] brs = new BufferedReader[2];
        brs[0] = IOUtils.getTextReader(vcfFile1);
        brs[1] = IOUtils.getTextReader(vcfFile2);
        try {
            String[] temps = new String[vcfFileNumber];
            String[][] tems = new String[vcfFileNumber][];
            for (int i = 0; i < vcfFileNumber; i++) {
                while ((temps[i] = brs[i].readLine()).startsWith("##")) {
                }
                tems[i] = temps[i].split("\t");
                for (int j = 0; j < tems[i].length; j++) {
                    int index = Collections.binarySearch(taxaListOlds[i], tems[i][j]);
                    if (index < 0) {
                        continue;
                    }
                    columnIndices[i][index] = j;
                }
            }
            int[] cnts = new int[2];
            for (int i = 0; i < allTaxaList.size(); i++) {
                int tIndex = columnIndices[tableIndex[i]][cnts[tableIndex[i]]];
                columnIndex[i] = tIndex;
                cnts[tableIndex[i]]++;
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(this.getVCFHeaderABD_AB(allTaxaList));
            bw.newLine();
            while ((temps[0] = brs[0].readLine()) != null) {
                temps[1] = brs[1].readLine();
                for (int i = 0; i < vcfFileNumber; i++) {
                    tems[i] = temps[i].split("\t");
                }
                StringBuilder sb = new StringBuilder();
                double maf = 0;
                double missing = 0;
                List<String> allGenoList = new ArrayList();
                List<String> abdGenoList = new ArrayList();
                List<String> abGenoList = new ArrayList();
                for (int i = 0; i < allTaxaList.size(); i++) {
                    if (tableIndex[i] == 0) {
                        abGenoList.add(tems[0][columnIndex[i]]);
                    } else {
                        abdGenoList.add(tems[1][columnIndex[i]]);
                    }
                    allGenoList.add(tems[tableIndex[i]][columnIndex[i]]);
                }
                String[] genoArray = allGenoList.toArray(new String[allGenoList.size()]);
                String[] hexaGenoArray = abdGenoList.toArray(new String[abdGenoList.size()]);
                String[] tetraGenoArray = abGenoList.toArray(new String[abGenoList.size()]);
                sb = new StringBuilder();
                for (int i = 0; i < 6; i++) {
                    sb.append(tems[0][i]).append("\t");
                }
                String info = this.getInfo(genoArray, tems[0][4]);
                sb.append(info);
                bw.write(sb.toString());
                bw.newLine();
                
                
                
                
                
            }
            
            
            for (int i = 0; i < vcfFileNumber; i++) {
                brs[i].close();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
    
    private String getInfo (String[] genoArray, String altList) {
        int dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1+nAlt];
        int[] acCnt = new int[1+nAlt];
        int[][] gnCnt = new int[1+nAlt][1+nAlt];
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":");
            temList = PStringUtils.fastSplit(tempList.get(1), ",");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                dp+=c;
                adCnt[j] += c;
            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                acCnt[c]++;
            }
            int index1 = Integer.parseInt(temList.get(0));
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++;
            if (index1 != index2) ht++;
        }
        nz = genoArray.length - nz;
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum+=acCnt[i];
        }
        double maf = ((double) acCnt[0] / sum);
        if (maf >= 0.5) {
            maf = ((double) acCnt[1] / sum);
        }
        StringBuilder sb = new StringBuilder();
        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            sb.append(adCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";AC=");
        for (int i = 1; i < acCnt.length; i++) {
            sb.append(acCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) {
            for (int j = i + 1; j < gnCnt.length; j++) {
                sb.append(gnCnt[i][j]).append(",");
            }
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";HT=").append(ht).append(";MAF=").append(maf);
        return sb.toString();
    }
    
    
    public String getVCFHeaderABD_AB(List<String> taxaList) {
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);
        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.1\n"
                + "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
        sb.append("##fileDate=").append(S.split(" ")[0]).append("\n");
        sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n"
                + "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n"
                + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
                + "##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n"
                + "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n"
                + "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n"
                + "##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n"
                + "##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n"
                + "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n"
                + "##INFO=<ID=AAF_ABD,Number=1,Type=Float,Description=\"Alternative allele frequency on hexaploid bread wheat\">\n"
                + "##INFO=<ID=AAF_AB,Number=1,Type=Float,Description=\"Alternative allele frequency on tetraploid emmer wheat\">\n"
                + "##ALT=<ID=D,Description=\"Deletion\">\n"
                + "##ALT=<ID=I,Description=\"Insertion\">\n"
                + "##Species=Wheat\n"
                + "##ReferenceGenome=iwgsc_refseqv1.0\n"
                + "##VariationMapVersion=vmap2\n");
        sb.append("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
        for (int i = 0; i < taxaList.size(); i++) {
            sb.append("\t").append(taxaList.get(i));
        }
        return sb.toString();
    }

    public void mergePosList() {
        String inFileS1 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergePosList/chr001.ABDgenome.10000lines.vcf.gz";
        String inFileS2 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergePosList/chr001.ABgenome.10000lines.vcf.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergePosList/out.txt";
        String[] alleles = {"A", "C", "G", "T", "D", "I"};
        Arrays.sort(alleles);
        double[] fre = new double[6];
        try {
            int chr1 = Integer.MIN_VALUE;
            int chr2 = Integer.MIN_VALUE;
            int taxaNum1 = Integer.MIN_VALUE;
            int taxaNum2 = Integer.MIN_VALUE;
            TIntArrayList posList1 = new TIntArrayList();
            List<String> referList1 = new ArrayList<>();
            List<String> altList1 = new ArrayList<>();
            List<String> altDepthList1 = new ArrayList<>();
            TIntArrayList posList2 = new TIntArrayList();
            List<String> referList2 = new ArrayList<>();
            List<String> altList2 = new ArrayList<>();
            List<String> altDepthList2 = new ArrayList<>();
            BufferedReader br = IOUtils.getTextGzipReader(inFileS1);
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")) {
            };
            taxaNum1 = temp.split("\t").length - 9;
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                temp = temp.substring(0, 150);
                tem = temp.split("\t");
                chr1 = Integer.parseInt(tem[0]);
                posList1.add(Integer.parseInt(tem[1]));
                referList1.add(tem[3]);
                altList1.add(tem[4]);
                altDepthList1.add(tem[7].split(";")[1].replace("AD=", ""));
            }
            br.close();
            br = IOUtils.getTextGzipReader(inFileS2);
            temp = null;
            while ((temp = br.readLine()).startsWith("##")) {
            };
            taxaNum2 = temp.split("\t").length - 9;
            double weight1 = (double) taxaNum1 / (taxaNum1 + taxaNum2);
            double weight2 = (double) taxaNum2 / (taxaNum1 + taxaNum2);
            tem = null;
            while ((temp = br.readLine()) != null) {
                temp = temp.substring(0, 50);
                tem = temp.split("\t");
                chr2 = Integer.parseInt(tem[0]);
                posList2.add(Integer.parseInt(tem[1]));
                referList2.add(tem[3]);
                altList2.add(tem[4]);
                altDepthList2.add(tem[7].split(";")[1].replace("AD=", ""));
            }
            if (chr1 != chr2) {
                System.out.println("Wrong input files! Program quits.");
                System.exit(0);
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tRef\tAlt");
            bw.newLine();
            TIntHashSet mergedPosSet = new TIntHashSet(posList1);
            mergedPosSet.addAll(posList2);
            int[] mergedPos = mergedPosSet.toArray();
            Arrays.sort(mergedPos);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < mergedPos.length; i++) {
                sb = new StringBuilder();
                int index1 = posList1.binarySearch(mergedPos[i]);
                int index2 = posList2.binarySearch(mergedPos[i]);
                if (index1 < 0 && index2 > -1) {
                    if (altList2.get(index2).length() > 3) {
                        continue;
                    }
                    sb.append(chr1).append("\t").append(posList2.get(index2)).append("\t").append(referList2.get(index2)).append("\t").append(altList2.get(index2));
                } else if (index1 > -1 && index2 < 0) {
                    if (altList1.get(index1).length() > 3) {
                        continue;
                    }
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t").append(altList1.get(index1));
                } else {
                    for (int j = 0; j < fre.length; j++) {
                        fre[j] = -1;
                    }
                    tem = altList1.get(index1).split(",");
                    String[] fretem = altDepthList1.get(index1).split(",");
                    double[] depth = new double[fretem.length];
                    double[] fre1 = new double[depth.length];
                    double sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        depth[j] = Integer.parseInt(fretem[j]);
                        sum += depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre1[j] = depth[j] / sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        fre[index] = fre1[j + 1] * weight1;
                    }

                    tem = altList2.get(index2).split(",");
                    fretem = altDepthList2.get(index2).split(",");
                    depth = new double[fretem.length];
                    double[] fre2 = new double[depth.length];
                    sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        depth[j] = Integer.parseInt(fretem[j]);
                        sum += depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre2[j] = depth[j] / sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        if (fre[index] < 0) {
                            fre[index] = fre2[j + 1] * weight2;
                        } else {
                            fre[index] = fre[index] + fre2[j + 1] * weight2;
                        }
                    }

                    int[] indices = PArrayUtils.getIndexByDescendingValue(fre);
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t");
                    for (int j = 0; j < 2; j++) {
                        if (fre[indices[j]] > 0) {
                            sb.append(alleles[indices[j]]).append(",");
                        } else {
                            break;
                        }
                    }
                    sb.deleteCharAt(sb.length() - 1);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
