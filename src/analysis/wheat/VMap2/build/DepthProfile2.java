package analysis.wheat.VMap2.build;

import pgl.AppUtils;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

class DepthProfile2 {

    public DepthProfile2 () {
//        this.mkTaxaBamDGenome();
//        this.mkPopDepInput();
//        this.mkPerlBatch();
        this.modifyPopDepAB();
    }

    public void modifyPopDepAB () {
        String sourceAB = "/Volumes/VMap2_Fei/popdep_vmap2/round_01/AB";
        String sourceABD = "/Volumes/VMap2_Fei/popdep_vmap2/round_01/ABD";
        String outAB = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/AB";
        String outABD = "/Volumes/VMap2_Fei/popdep_vmap2/round_02/ABD";
        this.modifyPopDep(sourceAB, outAB);
        this.modifyPopDep(sourceABD, outABD);
    }

    private void modifyPopDep (String sourceDirS, String outDirS) {
        new File(outDirS).mkdir();
        List<File> fList = IOUtils.getFileListInDirEndsWith(sourceDirS, ".gz");
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                List<String> l = null;
                String temp = null;
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(temp);
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    public void mkPerlBatch () {
//        java -Xmx100g -jar PlantGenetics.jar -a PopDep -p ./popdep_input/036_parameter_popdep.txt > log.txt
        String perlfile = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/run_popdep.pl";
        int[] chrIDs = RefV1Utils.getChrIDsOfSubgenomeD();
        StringBuilder sb = new StringBuilder();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(perlfile);
            for (int i = 0; i < chrIDs.length; i++) {
                sb.setLength(0);
                String currentChr = PStringUtils.getNDigitNumber(3, chrIDs[i]);
                String filename = currentChr + "_parameter_popdep.txt";
                String logname = currentChr + "_log.txt";
                sb.append("system(\"").append("java -Xmx100g -jar PlantGenetics.jar -a PopDep -p ./popdep_input/");
                sb.append(filename).append(" > ").append(logname).append("\");");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void mkPopDepInput () {
        String sampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/source/parameters_popdep.txt";
        String taxaBamFileS = "/data1/home/feilu/taxaBam.txt";
        String outDirS = "/data1/home/feilu/out";
        String samtoolsPath = "/data1/programs/samtools-1.8/samtools";
        String parametersDirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popdep_input";
        ArrayList<String> plinelist = new ArrayList<>();
        int[] chrIDs = RefV1Utils.getChrIDsOfSubgenomeD();
        for (int i = 0; i < chrIDs.length; i++) {
            plinelist.clear();
            plinelist.add(taxaBamFileS);
            plinelist.add(String.valueOf(chrIDs[i]));
            plinelist.add(String.valueOf(RefV1Utils.getChrIDLength(chrIDs[i])));
            plinelist.add(samtoolsPath);
            plinelist.add("32");
            String currentChr = PStringUtils.getNDigitNumber(3, chrIDs[i]);
            String resultFileS = new File (outDirS, currentChr+"_popdep.txt.gz").getAbsolutePath();
            String parameterFileS = new File (parametersDirS, currentChr+"_parameter_popdep.txt").getAbsolutePath();
            plinelist.add(resultFileS);
            AppUtils.creatParameterFile(sampleFileS, null, plinelist, parameterFileS);
        }

    }

    private void mkTaxaBamDGenome () {
        String dFileSYafei = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/source/yafei.txt";
        String dFileSXuebo = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/source/xuebo.txt";
        String dOutfileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth2/popdep_input/taxaBam.txt";
        RowTable<String> t = new RowTable<>(dFileSYafei);
        List<String> l = t.getColumn(1);
        HashSet<String> s = new HashSet<>(l);
        t = new RowTable<>(dFileSXuebo);
        l = t.getColumn(5);
        s.addAll(l);
        l = new ArrayList<>(s);
        Collections.sort(l);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(dOutfileS);
            bw.write("Taxa	Bams(A list of bams of the taxon, seperated by the delimiter of Tab)");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < l.size(); i++) {
                sb.setLength(0);
                sb.append("D").append(i+1).append("\t").append(l.get(i));
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
}
