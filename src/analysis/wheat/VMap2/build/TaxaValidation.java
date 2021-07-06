package analysis.wheat.VMap2.build;

import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class TaxaValidation {

    public TaxaValidation () {
//        this.profileSNPNumber();
        this.sampleSNPs();
    }

    public void sampleSNPs () {
        int sampleSize = 50000;
        String snpNumFileS = "/data1/home/feilu/VMap2_raw_snpCount.txt";
        RowTable<String> t = new RowTable(snpNumFileS);
        int[] chromSNPCounts = t.getColumnAsIntArray(1);
        String inDirS = "/data2/yafei/004_Vmap3/VCF/Raw_VCF/Vmap2_Out";
        String abOutfileS = "/data1/home/feilu/ab_sampleSite.vcf.gz";
        String dOutfileS = "/data1/home/feilu/d_sampleSite.vcf.gz";
        int[] abChrs = new int[] {1,2,3,4};
        int[] dChrs = new int[]{5,6};
        this.sampleSNPs(abChrs, inDirS, sampleSize, chromSNPCounts, abOutfileS);
        this.sampleSNPs(dChrs, inDirS, sampleSize, chromSNPCounts, dOutfileS);
    }

    private void sampleSNPs (int[] chrIDs, String vcfDirS, int sampleSize, int[] chromSNPCounts, String outfileS) {
        List<File> vcfList = IOUtils.getFileListInDirEndsWith(vcfDirS, ".gz");
        int[] snpIntervals  = new int[chrIDs.length];
        int sum = 0;
        for (int i = 0; i < chrIDs.length; i++) {
            sum+=chromSNPCounts[chrIDs[i]-1];
        }
        for (int i = 0; i < chrIDs.length; i++) {
            snpIntervals[i] = (int)((double)chromSNPCounts[chrIDs[i]-1]/((double)chromSNPCounts[chrIDs[i]-1]/sum*sampleSize));
        }
        List<String>[] resultLists = new ArrayList[chrIDs.length];
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < chrIDs.length; i++) {
            resultLists[i] = new ArrayList<>();
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            String currentChr = PStringUtils.getNDigitNumber(3, chrIDs[i]);
            String inVCF = null;

            for (int j = 0; j < vcfList.size(); j++) {
                if (vcfList.get(j).getName().contains(currentChr)) {
                    inVCF = vcfList.get(j).getAbsolutePath();
                    break;
                }
            }
            try {
                BufferedReader br = IOUtils.getTextGzipReader(inVCF);
                String temp = null;
                int cnt = -1;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        resultLists[i].add(temp);
                    }
                    else {
                        cnt++;
                        if (cnt%snpIntervals[i] != 0) continue;
                        resultLists[i].add(temp);
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println("Done with "+inVCF);
        });
        try {
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            for (int i = 0; i < chrIDs.length; i++) {
                for (int j = 0; j < resultLists[i].size(); j++) {
                    bw.write(resultLists[i].get(j));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void profileSNPNumber () {
        String inDirS = "/data2/yafei/004_Vmap3/VCF/Raw_VCF/Vmap2_Out";
        String outfileS = "/data1/home/feilu/VMap2_raw_snpCount.txt";
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        int[] counts = new int[fList.size()];
        fList.parallelStream().forEach(f -> {
            int currentChr = Integer.parseInt(f.getName().split("\\.")[0].replaceFirst("chr", ""));
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    cnt++;
                }
                br.close();
                counts[currentChr-1] = cnt;
                System.out.println("Done with "+ String.valueOf(currentChr));
            }
            catch (Exception e) {
                e.printStackTrace();
                System.out.println("Issue with "+ f.getName());
                System.exit(1);
            }
        });
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tSNPNumber");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < counts.length; i++) {
                sb.setLength(0);
                sb.append(i+1).append("\t").append(counts[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
