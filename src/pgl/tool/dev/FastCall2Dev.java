package pgl.tool.dev;



import pgl.app.fastCall2.FastCall2;
import pgl.app.fastCall2.IndividualCount;
import pgl.app.fastCall2.IndividualGenotype;
import pgl.app.fastCall2.VariationLibrary;

class FastCall2Dev {

    public FastCall2Dev () {
        this.withCommandLine();
//        this.miscellaneous();
    }

    public void miscellaneous () {
//        this.viewIndividualGenotype();
//        this.viewVariationLibrary();
//        this.viewIndiCounts();
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
//        this.variationDiscovery();
//        this.buildLibrary();
//        this.viewLibrary();
//        this.customizeLibrary();
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
        sb.append("-o ").append("/usr/local/bin/samtools ");
        String[] args = sb.toString().split(" ");
        new FastCall2(args);
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
        new FastCall2(args);
    }

    public void viewLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-mod ").append("vlib ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001.lib.gz ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/miscellaneous/1_1_200001.lib.txt ");
        String[] args = sb.toString().split(" ");
        new FastCall2(args);
    }

    public void customizeLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-mod ").append("clib ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001.lib.gz ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/custom_position.txt ");
        sb.append("-c ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001_sub.lib.gz ");
        String[] args = sb.toString().split(" ");
        new FastCall2(args);
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