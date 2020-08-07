package pgl.tool.dev;

public class PGLToolDevGo {

    public PGLToolDevGo () {
//        this.alignmentDev();
        this.genotypeDev();
//        this.popgenDev();
//        this.genomeAnnotationDev();
//        this.fastCallDev();
//        this.hapScannerDev();
//        this.speedCallDev();
//        this.popdepDev();
    }

    public void popdepDev () {
        new PopDepDev();
    }

    public void speedCallDev () {
        new SpeedCallDev();
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


}
