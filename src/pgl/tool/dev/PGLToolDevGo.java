package pgl.tool.dev;

public class PGLToolDevGo {

    public PGLToolDevGo () {
//        this.alignmentDev();
//        this.genotypeDev();
        this.fastCallDev();
//        this.speedCallDev();
    }

    public void speedCallDev () {
        new SpeedCallDev();
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
