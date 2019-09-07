/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PArrayUtils;

/**
 *
 * @author feilu
 */
public class VMapII {
    
    public VMapII () {
        this.mergePosList();
    }
    
    public void mergePosList () {
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
            while ((temp = br.readLine()).startsWith("##")) {};
            taxaNum1 = temp.split("\t").length-9;
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                temp = temp.substring(0, 50);
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
            while ((temp = br.readLine()).startsWith("##")) {};
            taxaNum2 = temp.split("\t").length-9;
            double weight1 = (double)taxaNum1/(taxaNum1+taxaNum2);
            double weight2 = (double)taxaNum2/(taxaNum1+taxaNum2);
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
                    if (altList2.get(index2).length() > 3) continue;
                    sb.append(chr1).append("\t").append(posList2.get(index2)).append("\t").append(referList2.get(index2)).append("\t").append(altList2.get(index2));
                }
                else if (index1 > -1 && index2 < 0) {
                    if (altList1.get(index1).length() > 3) continue;
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t").append(altList1.get(index1));
                }
                else  {
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
                        sum+=depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre1[j] = depth[j]/sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        fre[index] = fre1[j+1]*weight1;
                    }
                    
                    tem = altList2.get(index2).split(",");
                    fretem = altDepthList2.get(index2).split(",");
                    depth = new double[fretem.length];
                    double[] fre2 = new double[depth.length];
                    sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        depth[j] = Integer.parseInt(fretem[j]);
                        sum+=depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre2[j] = depth[j]/sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        if (fre[index] < 0) {
                            fre[index] = fre2[j+1]*weight2;
                        }
                        else {
                            fre[index] = fre[index]+fre2[j+1]*weight2;
                        }
                    }
                    
                    int[] indices = PArrayUtils.getIndexByDescendingValue(fre);
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t");
                    for (int j = 0; j < 2; j++) {
                        if (fre[indices[j]] > 0) {
                            sb.append(alleles[indices[j]]).append(",");
                        }
                        else {
                            break;
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
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
