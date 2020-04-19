package pgl.tool.dev;

public class PGLToolDevGo {

    public PGLToolDevGo () {
//        this.alignmentDev();
//        this.genotypeDev();
//        this.fastCallDev();
//        this.hapScannerDev();
//        this.speedCallDev();
        this.popdepDev();
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

    public void alignmentDev () {
        new AlignmentDev ();
    }

    public void genotypeDev () {
        new GenotypeDev();
    }

}
