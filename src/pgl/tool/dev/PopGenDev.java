package pgl.tool.dev;

import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Arrays;

class PopGenDev {

    public PopGenDev () {
//        this.mkSourceGenotype();
//        this.mkSourceAncestral();
        this.fd();
    }

    public void fd () {
//        this.fdRange();
        this.chromosomeFd();
    }

    public void chromosomeFd () {
        String genotypeFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popg/source/chr001_1_10000000_vmap2.1.vcf.gz";
        String ancestralFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popg/source/chr001_1_10000000_secer_hv_ancestral.txt.gz";
        GenotypeGrid gg = new GenotypeGrid(genotypeFileS, GenoIOFormat.VCF_GZ);
        for (int i = 0; i < gg.getTaxaNumber(); i++) {
            System.out.println(gg.getTaxonName(i));
        }
    }

    public void fdRange () {
        String genotypeFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popg/source/chr001_1_10000000_vmap2.1.vcf.gz";
        String ancestralFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popg/source/chr001_1_10000000_secer_hv_ancestral.txt.gz";
        GenotypeGrid gg = new GenotypeGrid(genotypeFileS, GenoIOFormat.VCF_GZ);
        TIntArrayList posList = new TIntArrayList();
        TByteArrayList ancestralList = new TByteArrayList();
        try {
            BufferedReader br = IOUtils.getTextGzipReader(ancestralFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                posList.add(Integer.parseInt(tem[1]));
                ancestralList.add(tem[2].getBytes()[0]);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int[] pos = posList.toArray();
        int siteNumber = gg.getSiteNumber();
        int cnt = 0;
        int cnt2 = 0;
        for (int i = 0; i < siteNumber; i++) {
            int currentPos = gg.getPosition(i);
            int index = Arrays.binarySearch(pos, currentPos);
            if (index < 0) continue;
            if ((byte)gg.getReferenceAlleleBase(i) == ancestralList.get(index) || (byte)gg.getAlternativeAlleleBase(i) == ancestralList.get(index)) cnt2++;
            cnt++;
        }
        System.out.println(siteNumber);
        System.out.println(cnt);
        System.out.println((double)cnt/siteNumber);
        System.out.println(cnt2);
        System.out.println((double)cnt2/siteNumber);
    }

    public void mkSourceAncestral () {
        String infileS = "/Volumes/VMap2_Fei/ancestral/chr001_secer_hv_ancestral.txt.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popg/source/chr001_1_10000000_secer_hv_ancestral.txt.gz";
        int regionStart = 1; //inclusive
        int regionEnd = 10_000_000; //inclusive
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = br.readLine();
            bw.write(temp);bw.newLine();
            String[] tem;
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                int current = Integer.parseInt(tem[1]);
                if (current < regionStart) continue;
                if (current > regionEnd) break;
                bw.write(temp);bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkSourceGenotype() {
        String genoInputFileS = "/Volumes/VMap2_Fei/vcf/004_vmap2.1/chr001_vmap2.1.vcf.gz";
        String genoOutputFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/popg/source/chr001_1_10000000_vmap2.1.vcf.gz";
        int regionStart = 1; //inclusive
        int regionEnd = 10_000_000; //inclusive
        try {
            BufferedReader br = IOUtils.getTextGzipReader(genoInputFileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(genoOutputFileS);
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")) {
                bw.write(temp);bw.newLine();
            }
            bw.write(temp);bw.newLine();
            String[] tem;
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t", 100);
                int current = Integer.parseInt(tem[1]);
                if (current < regionStart) continue;
                if (current > regionEnd) break;
                bw.write(temp);bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
