package analysis;

import analysis.envGWAS.EnvGWASGo;
import analysis.wheat.VMap1.WheatVMap1Go;
import analysis.wheat.VMap2.build.VMap2BuildGo;
import pgl.infra.utils.Benchmark;
import pgl.tool.dev.PGLToolDevGo;
import pgl.data.LabDataGo;
import analysis.wheat.VMap2.WheatVMap2Go;

public class PrivateEntrance {

    public PrivateEntrance (String[] args) {
//        this.analysis();
        this.data();
//        this.tool();
    }
    
    private void tool () {
        new PGLToolDevGo();
    }
    
    private void data () {
        new LabDataGo ();
    }
    
    private void analysis () {
        this.wheatVMap1Go();
//        this.wheatVMap2BuildGo();
//        this.wheatVMap2Go();
//        this.envGWASGo();
//        this.labPublicGo();
    }
    
    public void wheatVMap1Go () {
        new WheatVMap1Go();
    }

    public void wheatVMap2BuildGo () {
        new VMap2BuildGo();
    }

    public void wheatVMap2Go () {
        new WheatVMap2Go();
    }

    public void labPublicGo () {
        new LabDataGo();
    }

    public void envGWASGo() {
        new EnvGWASGo();
    }

    public static void main (String[] args) {
        long start = System.nanoTime();
        new PrivateEntrance(args);
        System.out.println("Program finished in " + Benchmark.getTimeSpanMinutes(start) + " minutes");
    }
}
