package pgl.tool.fix;

import java.io.BufferedReader;
import java.io.InputStreamReader;

import static pgl.infra.utils.Benchmark.getTimeSpanNanoseconds;

public class FastCall2Fix {

    public FastCall2Fix () {
        this.samtoolsBAQ();
    }

    public void samtoolsBAQ () {
        String command1 = "samtools mpileup -q 30 -Q 20 -f /Users/feilu/Documents/analysisL/softwareFix/pgl/fastCall2Fix/ref/chr036.fa /Users/feilu/Documents/analysisL/softwareFix/pgl/fastCall2Fix/bam/P01.chr036.bam -o /Users/feilu/Documents/analysisL/softwareFix/pgl/fastCall2Fix/mpileup/out.txt";
        String command2 = "samtools mpileup -B -q 30 -Q 20 -f /Users/feilu/Documents/analysisL/softwareFix/pgl/fastCall2Fix/ref/chr036.fa /Users/feilu/Documents/analysisL/softwareFix/pgl/fastCall2Fix/bam/P01.chr036.bam -o /Users/feilu/Documents/analysisL/softwareFix/pgl/fastCall2Fix/mpileup/out-B.txt";
        try {
            long start = System.nanoTime();
            Runtime rt = Runtime.getRuntime();
            Process p = rt.exec(command1);
            p.waitFor();
            long t1 = getTimeSpanNanoseconds(start);
            System.out.println("command1 takes "+ String.valueOf(t1));

            start = System.nanoTime();
            rt = Runtime.getRuntime();
            p = rt.exec(command2);
            p.waitFor();
            long t2 = getTimeSpanNanoseconds(start);
            System.out.println("command2 takes "+ String.valueOf(t2));
            System.out.println((double)t1/(double)t2);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
