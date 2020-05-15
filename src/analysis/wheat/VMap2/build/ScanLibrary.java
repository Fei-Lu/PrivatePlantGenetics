package analysis.wheat.VMap2.build;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;

public class ScanLibrary {

    public ScanLibrary () {
        this.mkLibrary();
    }

    public void mkLibrary () {
        String abVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/ab";
        String abdVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/abd";
        String dVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/d";

        String outDirS = "/Volumes/VMap2_Fei/variationLibrary";
        this.mkFromDirS(abVCFDirS, abdVCFDirS, outDirS);
        this.mkFromDirS(dVCFDirS, abdVCFDirS, outDirS);
    }

    public void mkFromDirS (String dirS1, String dirS2, String outDirS) {
        List<File> fList = IOUtils.getFileListInDir(dirS1);
        LongAdder counter = new LongAdder();
        fList.parallelStream().forEach(f -> {
            String inFileS2 = f.getName().split("\\.")[0]+".ABDgenome.vcf.gz";
            inFileS2 = new File(dirS2, inFileS2).getAbsolutePath();
            String outfileS = new File (outDirS, f.getName().split("\\.")[0]+".lib.txt.gz").getAbsolutePath();
            int snpCount = this.mkFrom2VCF(f.getAbsolutePath(), inFileS2, outfileS);
            counter.add(snpCount);
            System.out.println(f.getName());
        });
        System.out.println(dirS1);
        System.out.println(counter.longValue());
    }

    public int mkFrom2VCF (String inFileS1, String inFileS2, String outfileS) {
//        String inFileS1 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergePosList/chr001.ABDgenome.10000lines.vcf.gz";
//        String inFileS2 = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergePosList/chr001.ABgenome.10000lines.vcf.gz";
//        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/vmap2/mergePosList/out.txt";
        String[] alleles = {"A", "C", "G", "T", "D", "I"};
        Arrays.sort(alleles);
        double[] fre = new double[6];
        int totalCnt = 0;
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
                temp = temp.substring(0, 150);
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
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            bw.write("Chr\tPos\tRef\tAlt");
            bw.newLine();
            TIntHashSet mergedPosSet = new TIntHashSet(posList1);
            mergedPosSet.addAll(posList2);
            int[] mergedPos = mergedPosSet.toArray();
            totalCnt = mergedPos.length;
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

                    int[] indices = PArrayUtils.getIndicesByDescendingValue(fre);
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
        return totalCnt;
    }

}
