package analysis.wheat.VMap2.build;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import pgl.graph.r.DensityPlot;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.SiteVCF;
import pgl.infra.dna.genot.VCFUtils;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

class FilterSNPs {

    public FilterSNPs () {
//        this.profileFeatures();
//        this.filterSNPsStep1();

    }


    public void filterSNPsStep1 () {
        //only retain biallelic SNPs
        double hetThresh = 0.05;
        double nonmissingThresh = 0.8;
        int macThresh = 2;
        String taxon1 = "AB_071";
        String taxon2 = "AB_129";
        String inDirS = "/data2/yafei/004_Vmap3/VCF/Raw_VCF/Vmap2_Out";
        String outDirS = "/data1/home/feilu/filter1";
//        String inDirS = "/Users/feilu/Documents/analysisL/production/vmap2/taxaValidation/sampleGenotype";
//        String outDirS = "/Users/feilu/Documents/analysisL/production/vmap2/taxaValidation/filter1";
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS, ".gz");
        fList.parallelStream().forEach(f -> {
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            int removeIndex = -1;
            int addIndex = -1;
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                while ((temp = br.readLine()).startsWith("##")){}
                String[] tem = temp.split("\t");
                List<String> nameList = new ArrayList<>();
                for (int i = 9; i < tem.length; i++) {
                    if (tem[i].equals(taxon2))  {
                        removeIndex = i-9;
                        continue;
                    }
                    if (tem[i].equals(taxon1)) {
                        nameList.add("AB_071M");
                        addIndex = i-9;
                    }
                    else {
                        nameList.add(tem[i]);
                    }
                }
                bw.write(VCFUtils.getVCFAnnotation());
                bw.newLine();
                bw.write(VCFUtils.getVCFHeader(nameList.toArray(new String[nameList.size()])));
                bw.newLine();
                short[] depRemove;
                short[] depAdd;
                short[] depNew;
                while ((temp = br.readLine()) != null) {
                    SiteVCF sv = new SiteVCF(temp);
                    if (sv.getAlleleNumber() > 2) continue;
                    if (removeIndex > -1) {
                        depRemove = sv.getAlleleDepth(removeIndex);
                        depAdd = sv.getAlleleDepth(addIndex);
                        if (depRemove[0] < 0) {
                            for (int i = 0; i < depRemove.length; i++) {
                                depRemove[i] = 0;
                            }
                        }
                        if (depAdd[0] < 0) {
                            for (int i = 0; i < depAdd.length; i++) {
                                depAdd[i] = 0;
                            }
                        }
                        depNew = new short[depAdd.length];
                        for (int i = 0; i < depAdd.length; i++) {
                            depNew[i] = (short)(depAdd[i]+depRemove[i]);
                        }
                        int sum = 0;
                        int anyV = 0;
                        for (int i = 0; i < depNew.length; i++) {
                            if (depNew[i] < 0) anyV = -1;
                            sum+=depNew[i];
                        }
                        if (sum == 0 || anyV < 0) {
                            for (int i = 0; i < depNew.length; i++) {
                                depNew[i] = Short.MIN_VALUE;
                            }
                        }
                        sv.deleteTaxon(removeIndex);
                        sv.setAlleleDepth(addIndex, depNew);
                    }
                    if (sv.getTotalAllelePresence()[1] < macThresh) continue;
                    double ratio = (double)sv.getNonMissingTaxaNumber()/sv.getTaxaNumber();
                    if (ratio < nonmissingThresh) continue;
                    ratio = (double)sv.getHeterozygotesNumber()/sv.getTaxaNumber();
                    if (ratio > hetThresh) continue;
                    bw.write(sv.getVCFRecord());
                    bw.newLine();
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+ f.getName());
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    public void profileFeatures () {
        //based on the results, the nonmissing cutoff is set to >0.8, het set to <0.05
        String abVCF = "/Users/feilu/Documents/analysisL/production/vmap2/taxaValidation/sampleGenotype/ab_sampleSite.vcf.gz";
        String dVCF = "/Users/feilu/Documents/analysisL/production/vmap2/taxaValidation/sampleGenotype/d_sampleSite.vcf.gz";
        String outDirS = "/Users/feilu/Documents/analysisL/production/vmap2/filterSNPs/features_raw";
        this.profileFeatures(abVCF, outDirS, "AB");
        this.profileFeatures(dVCF, outDirS, "D");
    }

    private void profileFeatures (String infileS, String outdirS, String genomeType) {
        String hetOut = new File(outdirS, genomeType+"_heterozygosity.pdf").getAbsolutePath();
        String nonmissingOut = new File(outdirS, genomeType+"_nonMissing.pdf").getAbsolutePath();
        GenotypeGrid gg = new GenotypeGrid(infileS, GenoIOFormat.VCF_GZ);
        int snpNumber = gg.getSiteNumber();
        DoubleArrayList hetList = new DoubleArrayList();
        DoubleArrayList nonmissingList = new DoubleArrayList();
        for (int i = 0; i < snpNumber; i++) {
            hetList.add((double)gg.getHeterozygoteNumberBySite(i)/gg.getTaxaNumber());
            nonmissingList.add((double)gg.getNonMissingNumberBySite(i)/gg.getTaxaNumber());
        }
        double[] het = hetList.toDoubleArray();
        double[] nonmissing = nonmissingList.toDoubleArray();
        DensityPlot d = new DensityPlot(het);
        d.setSmoothN(50000);
        d.setTitle(genomeType+"_heterozygosity");
        d.setXLab("Heterozygosity");
        d.setYLab("Density");
        d.setXLim(0, 0.2);
        d.saveGraph(hetOut);
        d = new DensityPlot(nonmissing);
        d.setSmoothN(50000);
        d.setTitle(genomeType+"_nonMissing");
        d.setXLab("NonMissing_rate");
        d.setYLab("Density");
        d.setXLim(0, 1);
        d.saveGraph(nonmissingOut);
    }

}
