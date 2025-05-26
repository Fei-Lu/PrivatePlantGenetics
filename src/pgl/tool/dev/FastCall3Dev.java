package pgl.tool.dev;



import pgl.app.fastCall2.FastCall2;
import pgl.app.fastCall2.IndividualCount;
import pgl.app.fastCall2.IndividualGenotype;
import pgl.app.fastCall2.VariationLibrary;
import pgl.infra.utils.Benchmark;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import pgl.app.fastCall3.FastCall3;

class FastCall3Dev {

    public FastCall3Dev() {
        this.withCommandLine();
//        this.miscellaneous();
    }


    public void miscellaneous () {
//        this.viewIndividualGenotype();
//        this.viewVariationLibrary();
//        this.viewIndiCounts();
//        this.testPileup();
    }

    /*
    Samtools117 is finished in 62.42966 seconds
    Samtools121 is finished in 51.792007 seconds
    Bcf121 is finished in 0.01449925 seconds

    samtools 1.21 is 18% faster than 1.17
    Bcftools doesn't output the pileup format of each site. It is slow because it calculate much more than needed (counting alleles).
    The time is short here because the inputstream did not capture the output from command (due to "|")
     */
    public void testPileup () {
        String bamFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall3/pileupTest/bam/A_0001.bam";
        String referenceFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall3/pileupTest/ref/chr001.fa.gz";

//        String bamFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/bams/TW0060.sub.bam";
//        String referenceFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa";

        String samtools117Path = "/Users/feilu/Software/samtools-1.17/samtools";
        String samtools121Path = "/Users/feilu/Software/samtools-1.21/samtools";
        String bcf121Path = "/Users/feilu/Software/bcftools-1.21/bcftools";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall3/pileupTest/output/out.txt";
        int chrom = 1;
        String commandSamtools117 = null;
        StringBuilder sb = new StringBuilder();
        sb.append(samtools117Path).append(" mpileup --no-output-ends -B -q 30 -Q 20 ").append(" -f ").append(referenceFileS)
                .append(" ").append(bamFileS).append(" -r ").append(chrom);
        commandSamtools117 = sb.toString();
        System.out.println(sb.toString());
        String commandSamtools121 = null;
        sb.setLength(0);
        sb.append(samtools121Path).append(" mpileup --no-output-ends -B -q 30 -Q 20 ").append(" -f ").append(referenceFileS)
                .append(" ").append(bamFileS).append(" -r ").append(chrom);
        commandSamtools121 =sb.toString();
        System.out.println(sb.toString());

        String commandBcf121 = null;
        sb.setLength(0);

        sb.append("sh -c ").append(bcf121Path).append(" mpileup -B -q 30 -Q 20 ").append(" -f ").append(referenceFileS)
                .append(" ").append(bamFileS).append(" -r ").append(chrom).append(" | ").append(bcf121Path).append(" call -mv -Ov");
        commandBcf121 = sb.toString();
        System.out.println(sb.toString());
        try {
            Runtime rt = Runtime.getRuntime();
            long startTime = System.nanoTime();
            String temp = null;
            Process p = rt.exec(commandSamtools117);

            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((temp = br.readLine()) != null) {

            }
            br.close();
            p.waitFor();

            System.out.println("Samtools117 is finished in " + (float) Benchmark.getTimeSpanSeconds(startTime)+ " seconds");
            startTime = System.nanoTime();
            p = rt.exec(commandSamtools121);
            temp = null;
            br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((temp = br.readLine()) != null) {

            }
            br.close();
            p.waitFor();
            System.out.println("Samtools121 is finished in " + (float) Benchmark.getTimeSpanSeconds(startTime)+ " seconds");

            startTime = System.nanoTime();
            p = rt.exec(commandBcf121);
            temp = null;
            br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            br.close();
            p.waitFor();
            System.out.println("Bcf121 is finished in " + (float) Benchmark.getTimeSpanSeconds(startTime)+ " seconds");
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void viewIndiCounts () {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/gen/indiCounts/TW0060/1_1_200001.iac.gz";
        String oufileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/miscellaneous/1_1_200001.iac.txt";
        IndividualCount ic = new IndividualCount(infileS);
        ic.writeTextFile(oufileS);
    }

    public void viewVariationLibrary () {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001.lib.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/miscellaneous/1_1_200001.lib.txt";
        VariationLibrary vl = new VariationLibrary(infileS);
        vl.writeTextFileS(outfileS);
    }

    public void viewIndividualGenotype() {
        String inputfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/ing/TW0060/1_1_200001.ing.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/miscellaneous/TW0060_ing.txt";
        IndividualGenotype g = new IndividualGenotype (inputfileS);
        g.getAlleleChromPosition(1);
        g.writeTextFile(outfileS);
    }

    public void withCommandLine () {
        this.variationDiscovery();
        this.buildLibrary();
        this.viewLibrary();
        this.customizeLibrary();
        this.scanGenotype();
    }
    public void variationDiscovery() {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-mod ").append("disc ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/inputfile/taxaBamMap.txt ");
        sb.append("-c ").append("0 ");
        sb.append("-d ").append("30 ");
        sb.append("-e ").append("20 ");
        sb.append("-f ").append("2 ");
        sb.append("-g ").append("0.2 ");
        sb.append("-h ").append("3 ");
        sb.append("-i ").append("0.8 ");
        sb.append("-j ").append("0.35 ");
        sb.append("-k ").append("0.2 ");
        sb.append("-l ").append("1:1,200000 ");
        sb.append("-m ").append("32 ");
        sb.append("-n ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/ing ");
        sb.append("-o ").append("/Users/feilu/Software/samtools-1.21/samtools ");
        String[] args = sb.toString().split(" ");
        new FastCall3(args);
    }

    public void buildLibrary() {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-mod ").append("blib ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa ");
        sb.append("-b ").append("1:1,200000 ");
        sb.append("-c ").append("2 ");
        sb.append("-d ").append("32 ");
        sb.append("-e ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/ing ");
        sb.append("-f ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib ");
        String[] args = sb.toString().split(" ");
        new FastCall3(args);
    }

    public void viewLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-mod ").append("vlib ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001.lib.gz ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/miscellaneous/1_1_200001.lib.txt ");
        String[] args = sb.toString().split(" ");
        new FastCall3(args);
    }

    public void customizeLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-mod ").append("clib ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001.lib.gz ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/custom_position.txt ");
        sb.append("-c ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001_sub.lib.gz ");
        String[] args = sb.toString().split(" ");
        new FastCall3(args);
    }

    public void scanGenotype() {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-mod ").append("scan ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/inputfile/taxaBamMap.txt ");
        sb.append("-c ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001.lib.gz ");
        sb.append("-d ").append("1:1,200000 ");
        sb.append("-e ").append("0 ");
        sb.append("-f ").append("30 ");
        sb.append("-g ").append("20 ");
        sb.append("-h ").append("0.05 ");
        sb.append("-i ").append("/usr/local/bin/samtools ");
        sb.append("-j ").append("32 ");
        sb.append("-k ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/gen/ ");
        String[] args = sb.toString().split(" ");
        new FastCall2(args);
    }
}