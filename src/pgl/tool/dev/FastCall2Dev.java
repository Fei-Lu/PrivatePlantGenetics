package pgl.tool.dev;



import pgl.app.fastCall2.FastCall2;

class FastCall2Dev {

    public FastCall2Dev () {
//        this.callIndividual();
//        this.mkVariationLib();
        this.genotyping();
    }

    public void genotyping () {
        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/parameters_fastcall2_3.txt";
        new FastCall2(parameterFileS);
    }

    public void mkVariationLib () {
        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/parameters_fastcall2_2.txt";
        new FastCall2(parameterFileS);
    }

    public void callIndividual () {
        String parameterFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/fastCall2/inputfile/parameters_fastcall2_1.txt";
        new FastCall2(parameterFileS);
    }
}
