package pgl.tool.dev;

import pgl.app.fastCall.FastCall;
import pgl.infra.dna.*;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;

public class FastCallDev {

    public FastCallDev () {
//        this.mkSubReference();
//        this.mkBams();
//        this.mkTaxaBamFiles();
//        this.runFastCall();
    }

    public void runFastCall () {
        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/inputfile/parameters_fastcall.txt";
        new FastCall(parameterFileS);
    }

    public void mkTaxaBamFiles() {
        String bamDirS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/bams";
        String vmapIIFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/inputfile/taxaBamMap_VMapII.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/inputfile/taxaBamMap.txt";
        File[] fs = IOUtils.listFilesEndsWith(new File(bamDirS).listFiles(), ".bam");
        Arrays.sort(fs);
        String[] names = new String[fs.length];
        double[] depths = new double[fs.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = fs[i].getName().split("\\.")[0];
        }
        RowTable<String> t = new RowTable<>(vmapIIFileS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa	Coverage-Of-All-Bams	Bams(A list of bams of the taxon, seperated by the delimiter of Tab)");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(names, t.getCell(i, 0));
                if (index < 0) continue;
                sb.setLength(0);
                sb.append(t.getCell(i, 0)).append("\t").append(t.getCell(i,1)).append("\t").append(fs[index].getAbsolutePath());
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

    public void mkBams () {
        //samtools view -b /data3/wgs/bam/ABD/TW0060.rmdup.bam 1:1-1000000 -o TW0060.sub.bam
        String indirS = "/data3/wgs/bam/ABD/";
        int start = 60;
        int end = 100;
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/bams/bamMk.pl";
        String outfileS2 = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/bams/bamIndex.pl";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 60; i < 100; i++) {
                sb.setLength(0);
                sb.append("system(\"samtools view -b ").append(indirS).append("TW00").append(i).append(".rmdup.bam 1:1-1000000 -o TW00").append(i).append(".sub.bam &\");");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IOUtils.getTextWriter(outfileS2);
            for (int i = 60; i < 100; i++) {
                sb.setLength(0);
                sb.append("system(\"samtools index TW00").append(i).append(".sub.bam &\");");
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

    public void mkSubReference () {
        int subGenomeSize = 1_000_000;
        String originalRef = "/Users/feilu/Documents/database/wheat/reference/v1.0/byChr/chr001.fa.gz";
        String targetRef = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa.gz";
        FastaBit fa = new FastaBit(originalRef);
        FastaRecordBit fb = fa.getFastaRecordBit(0, 0, subGenomeSize);
        fa = new FastaBit(fb);
        fa.writeFasta(targetRef, IOFileFormat.TextGzip);
    }

}

