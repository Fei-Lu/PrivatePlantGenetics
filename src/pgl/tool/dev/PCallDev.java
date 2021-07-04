package pgl.tool.dev;

import pgl.app.fastCall2.FastCall2;
import pgl.graph.r.Histogram;
import pgl.infra.table.RowTable;

public class PCallDev {

    public PCallDev () {
//        this.callIndividual();
//        this.mkVariationLib();
//        this.genotyping();
        depthTest();
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

    public void depthTest () {
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/pCall/depth/a.txt";
        RowTable<String> t = new RowTable<>(infileS);
        double[] values = t.getColumnAsDoubleArray(2);
        Histogram h  = new Histogram(values);
        h.showGraph();
    }
}
