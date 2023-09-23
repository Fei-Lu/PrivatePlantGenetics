package pgl.tool.dev;

import pgl.app.hapScanner.HapScanner;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeRows;
import pgl.infra.dna.genot.GenotypeTable;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;

public class HapScannerDev {

    public HapScannerDev () {
//        this.mkTaxaRefBamFileS();
//        this.mkPosAlleleFileS();
//        this.mkPosFileS();
        this.runHapScanner();
    }

    public void runHapScanner () {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/hapScanner/parameters_hapScanner.txt";
        new HapScanner(infileS);
    }

    public void mkPosFileS() {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/hapScanner/inputfile/posAllele_hapscanner.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/hapScanner/inputfile/pos_hapscanner.txt";
        RowTable<String> t = new RowTable(infileS);
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append(t.getCell(i,0)).append("\t").append(t.getCell(i,1));
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

    public void mkPosAlleleFileS () {
        String vcfFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/vcf/chr001.vcf";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/hapScanner/inputfile/posAllele_hapscanner.txt";
        int size = 50;
        GenotypeTable gt = new GenotypeRows(vcfFileS, GenoIOFormat.VCF);
        StringBuilder sb = new StringBuilder();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tRef\tAlt(maximum 2 alternative alleles, which is seperated by \",\", e.g. A,C)");
            bw.newLine();
            for (int i = 0; i < size; i++) {
                sb.setLength(0);
                sb.append(gt.getChromosome(i)).append("\t").append(gt.getPosition(i)).append("\t").append(AlleleEncoder.getAlleleBaseFromCoding(gt.getReferenceAlleleByte(i)));
                sb.append("\t").append(AlleleEncoder.getAlleleBaseFromCoding(gt.getAlternativeAlleleByte(i)));
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

    public void mkTaxaRefBamFileS() {
        String bamDirS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/bams";
        String ref = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/hapScanner/inputfile/taxaRefBAM_hapscanner.txt";
        int size = 10;
        File[] fs = new File(bamDirS).listFiles();
        Arrays.sort(fs);
        if (fs.length < size) size = fs.length;
        fs = IOUtils.listFilesEndsWith(fs, ".bam");
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tReference\tBamPath");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            String name = null;
            for (int i = 0; i < size; i++) {
                sb.setLength(0);
                name = fs[i].getName().substring(0,6);
                sb.append(name).append("\t").append(ref).append("\t").append(fs[i].getAbsolutePath());
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
