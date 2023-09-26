package entrance;

import analysis.AnalysisSection;
import learning.LearnSection;
import pgl.infra.utils.Benchmark;
import pgl.tool.dev.PGLToolSection;
import pgl.data.PGLDataSection;

public class PrivateEntrance {

    public PrivateEntrance (String[] args) {
//        this.analysis();
//        this.data();
//        this.tool();
        this.learning();
    }

    private void learning () {
        new LearnSection();
    }

    private void tool () {
        new PGLToolSection();
    }

    private void data() {
        new PGLDataSection();
    }

    private void analysis () {
        new AnalysisSection();
    }

    public static void main (String[] args) {
        long start = System.nanoTime();
        new PrivateEntrance(args);
        System.out.println("Program is finished in " + (float)Benchmark.getTimeSpanMinutes(start) + " minutes");
    }
}
