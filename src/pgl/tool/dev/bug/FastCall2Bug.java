package pgl.tool.dev.bug;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.app.fastCall2.VariationLibrary;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.sql.PreparedStatement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

class FastCall2Bug {

    public FastCall2Bug () {
        this.nAltAllele();
    }

    public void nAltAllele () {
//        this.findNAlt();
//        this.checkOriginalVMap2();
//        this.compareN();
        this.checkLibrary();
    }

    private void checkLibrary () {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/1_1_471304006.lib.gz";
        String positionFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/nAlt.beforeAfter.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/1_1_471304006.lib.txt";
        VariationLibrary vl = new VariationLibrary(infileS);
        RowTable<String> t = new RowTable<>(positionFileS);
        int[] poses = t.getColumnAsIntArray(1);
        int[] positionIndices = new int[poses.length];
        for (int i = 0; i < poses.length; i++) {
            positionIndices[i] = vl.getPositionIndex(poses[i]);
        }
        vl.writeTextFileS(outfileS, positionIndices);
    }

    private void compareN () {
        String infileS1 = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/nAlt.info.txt";
        String infileS2 = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/nAlt.oVMap2.info.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/nAlt.beforeAfter.txt";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS2);
            IntArrayList posList = new IntArrayList();
            List<String> altList = new ArrayList<>();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                posList.add(Integer.parseInt(tem[1]));
                altList.add(tem[4]);
            }
            br.close();
            br = IOUtils.getTextReader(infileS1);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tRef\tAltBefore\tAltAfter");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            int[] poses = posList.toIntArray();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                int currentPos = Integer.parseInt(tem[1]);
                int index = Arrays.binarySearch(poses, currentPos);
                if (index < 0) continue;
                sb.setLength(0);
                sb.append(tem[0]).append("\t").append(tem[1]).append("\t")
                        .append(tem[3]).append("\t").append(altList.get(index)).append("\t").append(tem[4]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void checkOriginalVMap2 () {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/nAlt.info.txt.gz";
        String vcfFileS = "/Volumes/VMap2_Fei/chr001_vmap2.0.vcf.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/nAlt.oVMap2.info.txt";
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            IntArrayList posList = new IntArrayList();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                posList.add(Integer.parseInt(l.get(1)));
            }
            br.close();
            int[] poses = posList.toIntArray();
            br = IOUtils.getTextGzipReader(vcfFileS);
            while ((temp = br.readLine()).startsWith("##")) continue;
            int cntN = 0;
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt));
                l = PStringUtils.fastSplit(temp.substring(0, 200));
                int currentPos = Integer.parseInt(l.get(1));
                if (Arrays.binarySearch(poses, currentPos) < 0) continue;
                sb.setLength(0);
                for (int i = 0; i < 9; i++) {
                    sb.append(l.get(i)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
                cntN++;
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(String.valueOf(cntN));
            System.out.println(String.valueOf(cnt));
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void findNAlt () {
        String infileS = "/Volumes/VMap2_Fei/vcf/round_02/001_raw_vcf/chr001.vmap2.vcf.gz";
        String outFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/bugFix/nAlt/nAlt.info.txt";
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()).startsWith("##")) continue;
            List<String> l = new ArrayList<>();
            StringBuilder sb = new StringBuilder();
            int cntN = 0;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt));
                l = PStringUtils.fastSplit(temp.substring(0, 200));
                if (!l.get(4).contains("N")) continue;
                sb.setLength(0);
                for (int i = 0; i < 9; i++) {
                    sb.append(l.get(i)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
                cntN++;
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(String.valueOf(cntN));
            System.out.println(String.valueOf(cnt));
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
