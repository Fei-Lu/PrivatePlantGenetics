package pgl.tool.dev;



import pgl.app.fastCall2.FastCall2;

class FastCall2Dev {

    public FastCall2Dev () {
//        this.withParameterFileS();
        this.withCommandLine();
    }

    public void test () {

    }

    public void withCommandLine () {
        this.commandDiscovery();
        this.commandLibrary();
        this.commandGenotype();
    }
    public void commandDiscovery () {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-tool ").append("disc ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/inputfile/taxaBamMap.txt ");
        sb.append("-c ").append("30 ");
        sb.append("-d ").append("20 ");
        sb.append("-e ").append("2 ");
        sb.append("-f ").append("0.2 ");
        sb.append("-g ").append("3 ");
        sb.append("-h ").append("0.8 ");
        sb.append("-i ").append("0.35 ");
        sb.append("-j ").append("0.2 ");
        sb.append("-k ").append("1:1,200000 ");
        sb.append("-l ").append("32 ");
        sb.append("-m ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/ing ");
        sb.append("-n ").append("/usr/local/bin/samtools ");
        String[] args = sb.toString().split(" ");
        new FastCall2(args);
    }

    public void commandLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-tool ").append("blib ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa ");
        sb.append("-b ").append("1:1,200000 ");
        sb.append("-c ").append("2 ");
        sb.append("-d ").append("32 ");
        sb.append("-e ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/ing ");
        sb.append("-f ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib ");
        String[] args = sb.toString().split(" ");
        new FastCall2(args);
    }

    public void commandGenotype () {
        StringBuilder sb = new StringBuilder();
        sb.append("-app ").append("FastCall2 ");
        sb.append("-tool ").append("scan ");
        sb.append("-a ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/ref/chr001_1Mb.fa ");
        sb.append("-b ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall/inputfile/taxaBamMap.txt ");
        sb.append("-c ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/vLib/1_1_200001.lib.gz ");
        sb.append("-d ").append("1:1,200000 ");
        sb.append("-e ").append("0.05 ");
        sb.append("-f ").append("/usr/local/bin/samtools ");
        sb.append("-g ").append("32 ");
        sb.append("-h ").append("/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/gen/ ");
        String[] args = sb.toString().split(" ");
        new FastCall2(args);
    }

//    /**
//     * @deprecated
//     */
//    public void withParameterFileS () {
//        this.callIndividual();
//        this.mkVariationLib();
//        this.genotyping();
//    }
//
//    public void genotyping () {
//        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/parameters_fastcall2_3.txt";
//        new FastCall2(parameterFileS);
//    }
//
//    public void mkVariationLib () {
//        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/parameters_fastcall2_2.txt";
//        new FastCall2(parameterFileS);
//    }
//
//    public void callIndividual () {
//        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/parameters_fastcall2_1.txt";
//        new FastCall2(parameterFileS);
//    }
}
