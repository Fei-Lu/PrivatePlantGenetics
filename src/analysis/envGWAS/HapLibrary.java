package analysis.envGWAS;

import pgl.infra.utils.IOUtils;

import java.io.File;
import java.util.List;

public class HapLibrary {

    public HapLibrary () {
        this.pileupDepth();
    }

    public void pileupDepth () {
        String bamDirS = "/Users/feilu/Documents/analysisH/envGWAS/bam";
        String samtools = "/usr/local/bin/samtools";
        String refFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/ABD/abd_iwgscV1.fa";
        String posFileS = "/Users/feilu/Documents/analysisH/envGWAS/pos.txt";
        String outFileS = "/Users/feilu/Documents/analysisH/envGWAS/bam.pileup.txt";
        List<File> fList = IOUtils.getFileListInDirEndsWith(bamDirS, "bam");
        StringBuilder sb = new StringBuilder(samtools);
        sb.append(" mpileup -A -B -q 20 -Q 20 -f ").append(refFileS).append(" ");
        for (int i = 0; i < fList.size(); i++) {
            sb.append(fList.get(i).getAbsolutePath()).append(" ");
        }
        sb.append("-l ").append(posFileS).append(" -o ").append(outFileS);
        System.out.println(sb.toString());
        fList.parallelStream().forEach(f -> {
            System.out.println(f.getAbsoluteFile());
        });
//        try {
//            Runtime rt = Runtime.getRuntime();
//            Process p = rt.exec(sb.toString());
//            p.waitFor();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
        System.out.println("Done pileup");
    }
}
