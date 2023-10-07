package pgl.tool.dev;

import pgl.tool.fix.FastCall2Fix;

public class PGLToolSection {

    public PGLToolSection() {
        this.dev();
//        this.fix();
    }

    public void fix() {
        new FastCall2Fix();
    }

    public void dev () {
//        this.dnaDev();
//        this.alignmentDev();
//        this.genotypeDev();
//        this.popgenDev();
//        this.genomeAnnotationDev();
//        this.fastCallDev();
//        this.hapScannerDev();
//        this.fastCall2Dev();
        this.graphicsDev();
//        this.popdepDev();
    }

    public void popdepDev () {
        new PopDepDev();
    }

    public void graphicsDev() {
        new GraphicsDev();
    }

    public void fastCall2Dev () {
        new FastCall2Dev();
    }

    public void hapScannerDev () {
        new HapScannerDev();
    }

    public void fastCallDev () {
        new FastCallDev();
    }

    public void genomeAnnotationDev () {
        new GenomeAnnotationDev();
    }

    public void popgenDev () {
        new PopGenDev();
    }

    public void alignmentDev () {
        new AlignmentDev ();
    }

    public void genotypeDev () {
        new GenotypeDev();
    }

    public void dnaDev () {
        new DNADev();
    }

}
