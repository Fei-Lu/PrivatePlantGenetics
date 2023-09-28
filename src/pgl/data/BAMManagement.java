package pgl.data;

import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;

public class BAMManagement {

    public BAMManagement () {
        //this.listAllBams();
//        this.mkNewpaths();
//        this.mkMVCommand();
    }

    public void mkMVCommand () {
        String infileS =  "/Users/feilu/Documents/analysisL/production/bamManage/bam_Old2New.map.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/production/bamManage/moveBam.sh";
        RowTable<String> t = new RowTable<>(infileS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                sb.append("mv ").append(t.getCell(i, 0)).append(" ").append(t.getCell(i, 1));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mkNewpaths () {
        String infileS = "/Users/feilu/Documents/analysisL/production/bamManage/oldPaths.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/production/bamManage/bam_Old2New.map.txt";
        String aDirS = "/data3/wgs/bam/A/";
        String abDirS = "/data3/wgs/bam/AB/";
        String abdDirS = "/data3/wgs/bam/ABD/";
        String dDirS = "/data3/wgs/bam/D/";
        RowTable<String> t = new RowTable<>(infileS);
        int aCountBam = 1;
        int abCountBam = 1;
        int abdCountBam = 1;
        int dCountBam = 1;
        int aCountBai = 1;
        int abCountBai = 1;
        int abdCountBai = 1;
        int dCountBai = 1;
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("OldPath\tNewPath");
            bw.newLine();
            String newName = null;
            StringBuilder sb = new StringBuilder();;
            for (int i = 0; i < t.getRowNumber(); i++) {
                String parent = new File (t.getCell(i,0)).getParent();
                String oldName = new File (t.getCell(i,0)).getName();
                sb.setLength(0);
                if (parent.endsWith("/A")) {
                    if (oldName.endsWith(".bam")) {
                        newName = "A_"+ PStringUtils.getNDigitNumber(4, aCountBam)+".bam";
                        aCountBam++;
                    }
                    else if (oldName.endsWith(".bai")) {
                        newName = "A_"+ PStringUtils.getNDigitNumber(4, aCountBai)+".bam.bai";
                        aCountBai++;
                    }
                    sb.append(t.getCell(i,0)).append("\t").append(new File(aDirS, newName).getAbsolutePath());
                }
                else if (parent.endsWith("/AB")) {
                    if (oldName.endsWith(".bam")) {
                        newName = "AB_"+ PStringUtils.getNDigitNumber(4, abCountBam)+".bam";
                        abCountBam++;
                    }
                    else if (oldName.endsWith(".bai")) {
                        newName = "AB_"+ PStringUtils.getNDigitNumber(4, abCountBai)+".bam.bai";
                        abCountBai++;
                    }
                    sb.append(t.getCell(i,0)).append("\t").append(new File(abDirS, newName).getAbsolutePath());
                }
                else if (parent.endsWith("/ABD")) {
                    if (oldName.startsWith("CS")) {
                        if (oldName.endsWith(".bam")) {
                            newName = "CS_sg_2014_3X"+ ".bam";
                        }
                        else if (oldName.endsWith(".bai")) {
                            newName = "CS_sg_2014_3X"+".bam.bai";
                        }
                        sb.append(t.getCell(i,0)).append("\t").append(new File(abdDirS, newName).getAbsolutePath());
                    }
                    else {
                        if (oldName.endsWith(".bam")) {
                            newName = "ABD_"+ PStringUtils.getNDigitNumber(4, abdCountBam)+".bam";
                            abdCountBam++;
                        }
                        else if (oldName.endsWith(".bai")) {
                            newName = "ABD_"+ PStringUtils.getNDigitNumber(4, abdCountBai)+".bam.bai";
                            abdCountBai++;
                        }
                        sb.append(t.getCell(i,0)).append("\t").append(new File(abdDirS, newName).getAbsolutePath());
                    }
                }
                else if (parent.endsWith("/D")) {
                    if (oldName.endsWith(".bam")) {
                        newName = "D_"+ PStringUtils.getNDigitNumber(4, dCountBam)+".bam";
                        dCountBam++;
                    }
                    else if (oldName.endsWith(".bai")) {
                        newName = "D_"+ PStringUtils.getNDigitNumber(4, dCountBai)+".bam.bai";
                        dCountBai++;
                    }
                    sb.append(t.getCell(i,0)).append("\t").append(new File(dDirS, newName).getAbsolutePath());
                }
                else if (parent.endsWith("rmdupbam")){
                    if (oldName.startsWith("CS")) {
                        if (oldName.endsWith(".bam")) {
                            newName = "CS_mp_2018_8X"+ ".bam";
                        }
                        else if (oldName.endsWith(".bai")) {
                            newName = "CS_mp_2018_8X"+".bam.bai";
                        }
                        sb.append(t.getCell(i,0)).append("\t").append(new File(abdDirS, newName).getAbsolutePath());
                    }
                }
                else if (parent.endsWith("ab_S25")){//AB genome
                    if (oldName.endsWith(".bam")) {
                        newName = "AB_"+ PStringUtils.getNDigitNumber(4, abCountBam)+".bam";
                        abCountBam++;
                    }
                    else if (oldName.endsWith(".bai")) {
                        newName = "AB_"+ PStringUtils.getNDigitNumber(4, abCountBai)+".bam.bai";
                        abCountBai++;
                    }
                    sb.append(t.getCell(i,0)).append("\t").append(new File(abDirS, newName).getAbsolutePath());
                }
                else if (parent.endsWith("d_S5")){//D genome
                    if (oldName.endsWith(".bam")) {
                        newName = "D_"+ PStringUtils.getNDigitNumber(4, dCountBam)+".bam";
                        dCountBam++;
                    }
                    else if (oldName.endsWith(".bai")) {
                        newName = "D_"+ PStringUtils.getNDigitNumber(4, dCountBai)+".bam.bai";
                        dCountBai++;
                    }
                    sb.append(t.getCell(i,0)).append("\t").append(new File(dDirS, newName).getAbsolutePath());
                }
                else if (parent.contains("abd_S61")){
                    if (oldName.endsWith(".bam")) {
                        newName = "ABD_"+ PStringUtils.getNDigitNumber(4, abdCountBam)+".bam";
                        abdCountBam++;
                    }
                    else if (oldName.endsWith(".bai")) {
                        newName = "ABD_"+ PStringUtils.getNDigitNumber(4, abdCountBai)+".bam.bai";
                        abdCountBai++;
                    }
                    sb.append(t.getCell(i,0)).append("\t").append(new File(abdDirS, newName).getAbsolutePath());
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
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
                String name = fs[i].getName();
                if (name.endsWith(".bam") || name.endsWith(".bai")) {
                    bw.write(fs[i].getAbsolutePath());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
