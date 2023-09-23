package analysis;

import analysis.envGWAS.EnvGWASGo;
import analysis.fastcall.FastCallGo;
import analysis.wheat.RNASeq.WEGAGo;
import analysis.wheat.VMap1.WheatVMap1Go;
import analysis.wheat.VMap2.WheatVMap2Go;
import analysis.wheat.VMap2.build.VMap2BuildGo;
import analysis.wheat.VMap3.WheatVMap3Go;
import analysis.wheat.sv.WheatSVGo;

public class AnalysisSection {
    public AnalysisSection () {
        this.wegaGo();
        this.svGo();
        this.wheatVMap1Go();
        this.wheatVMap2BuildGo();
        this.wheatVMap2Go();
        this.wheatVMap3Go();
        this.envGWASGo();
        this.fastCallGo();
    }


    public void fastCallGo () {
        new FastCallGo();
    }

    public void wegaGo () {
        new WEGAGo();
    }

    public void svGo () {
        new WheatSVGo();
    }

    public void wheatVMap1Go () {
        new WheatVMap1Go();
    }

    public void wheatVMap2BuildGo () {
        new VMap2BuildGo();
    }

    public void wheatVMap3Go () {
        new WheatVMap3Go();
    }

    public void wheatVMap2Go () {
        new WheatVMap2Go();
    }

    public void envGWASGo() {
        new EnvGWASGo();
    }

}
