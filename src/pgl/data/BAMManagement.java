package pgl.data;

import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;

public class BAMManagement {

    public BAMManagement () {
        this.listAllBams();
    }


    public void listAllBams () {
        String infileDirS = "/data3/wgs/bam";
        String outfileS = "./oldPaths.txt";
        File[] fs = IOUtils.listRecursiveFiles(new File(infileDirS));
        Arrays.sort(fs);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("OldPaths");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                bw.write(fs[i].getAbsolutePath());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
