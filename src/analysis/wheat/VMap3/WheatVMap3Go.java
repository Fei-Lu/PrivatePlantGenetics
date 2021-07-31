package analysis.wheat.VMap3;

import pgl.infra.dna.genot.SiteVCF;
import pgl.infra.dna.genot.VCFUtils;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.StreamCorruptedException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class WheatVMap3Go {

    public WheatVMap3Go () {
        this.replaceSamples();
    }

    public void replaceSamples () {
        String inDirS1 = "";
        String inDirS2 = "";
        String outDirS = "";
        int[] chrIDsA = RefV1Utils.getChrIDsOfSubgenomeA();
        int[] chrIDsB = RefV1Utils.getChrIDsOfSubgenomeB();
        int[] chrIDs = new int[chrIDsA.length+chrIDsB.length];
        System.arraycopy(chrIDsA,0,chrIDs, 0, chrIDsA.length);
        System.arraycopy(chrIDsB, 0, chrIDs, chrIDsA.length, chrIDsB.length);
        Arrays.sort(chrIDs);
        List<Integer> chrIDList = new ArrayList<>();
        for (int i = 0; i < chrIDs.length; i++) {
            chrIDList.add(chrIDs[i]);
        }
        chrIDList.parallelStream().forEach(chr -> {
            String file1S = new File(inDirS1,PStringUtils.getNDigitNumber(3, chr)+"abc.vcf.gz").getAbsolutePath();
            String file2S = new File(inDirS2,PStringUtils.getNDigitNumber(3, chr)+"abc.vcf.gz").getAbsolutePath();
            String outfileS = new File(outDirS,PStringUtils.getNDigitNumber(3, chr)+"abc.vcf.gz").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(file1S);
                BufferedReader br2 = IOUtils.getTextGzipReader(file2S);
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp1 = null;
                String temp2 = null;
                while ((temp1 = br1.readLine()).startsWith("##")) {}
                while ((temp2 = br2.readLine()).startsWith("##")) {}
                String[] taxa1 = VCFUtils.getTaxaNames(temp1);
                String[] taxa2 = VCFUtils.getTaxaNames(temp2);
                int[] mapIndices = new int[taxa2.length];
                for (int i = 0; i < taxa2.length; i++) {
                    mapIndices[i] = Arrays.binarySearch(taxa1, taxa2[i]);
                }
                bw.write(VCFUtils.getVCFAnnotation());
                bw.newLine();
                bw.write(VCFUtils.getVCFHeader(taxa1));
                bw.newLine();
                SiteVCF vcf1 = null;
                SiteVCF vcf2 = null;
                int cnt = 0;
                while ((temp1 = br2.readLine()) != null) {
                    temp2 = br1.readLine();
                    vcf1 = new SiteVCF(temp1);
                    vcf2 = new SiteVCF(temp2);
                    for (int i = 0; i < taxa2.length; i++) {
                        vcf1.setAlleleDepth(mapIndices[i], vcf2.getAlleleDepth(i));
                    }
                    bw.write(vcf1.getVCFRecord());
                    bw.newLine();
                    cnt++;
                    if (cnt%2000000 == 0) System.out.println(String.valueOf(cnt)+" "+ outfileS);
                }
                bw.flush();
                bw.close();
                br1.close();
                br2.close();
                System.out.println(outfileS);
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }
}
