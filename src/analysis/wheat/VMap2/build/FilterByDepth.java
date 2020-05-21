package analysis.wheat.VMap2.build;

import gnu.trove.list.array.TDoubleArrayList;
import pgl.graphcis.r.Histogram;
import pgl.infra.dna.genotype.*;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;

class FilterByDepth {

    public FilterByDepth () {
//        this.filterVCF();
//        this.qualityCheck();
    }

    public void qualityCheck () {
        //missing, maf, heterozygous proportion
        String qcDirS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/qualityCheck";

        String abInVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/ab";
        String abdInVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/abd";
        String dInVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/d";

        this.checkQuality(abInVCFDirS, new File(qcDirS, "ab").getAbsolutePath(), "AB");
        this.checkQuality(abdInVCFDirS, new File(qcDirS, "abd").getAbsolutePath(), "ABD");
        this.checkQuality(dInVCFDirS, new File(qcDirS, "d").getAbsolutePath(), "D");

    }

    private void checkQuality (String vcfDirS, String outDirS, String genomeType) {
        File outDir = new File (outDirS);
        outDir.mkdir();
        int fN = 2;
        int size = 10000;
        List<File> fList = IOUtils.getFileListInDir(vcfDirS);
        Collections.sort(fList);

        GenotypeGrid[] gts = new GenotypeGrid[fN];
        int[] siteCounts = new int[gts.length];
        int totalSiteCount = 0;
        for (int i = 0; i < gts.length; i++) {
            gts[i] = new GenotypeGrid(fList.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
            totalSiteCount+=gts[i].getSiteNumber();
        }
        GenotypeOperation.mergeGenotypesBySite(gts[0], gts[1]);
        GenotypeGrid gt = gts[0];
        TDoubleArrayList missingSite = new TDoubleArrayList();
        TDoubleArrayList hetSite = new TDoubleArrayList();
        TDoubleArrayList missingTaxon = new TDoubleArrayList();
        TDoubleArrayList hetTaxon = new TDoubleArrayList();
        TDoubleArrayList maf = new TDoubleArrayList();
        int step = gt.getSiteNumber()/size;
        for (int i = 0; i < gt.getSiteNumber(); i+=step) {
            missingSite.add(((double) gt.getMissingNumberBySite(i)/gt.getTaxaNumber()));
            hetSite.add(gt.getHeterozygousProportionBySite(i));
            maf.add(gt.getMinorAlleleFrequency(i));
        }
        for (int i = 0; i < gt.getTaxaNumber(); i++) {
            missingTaxon.add((double)gt.getMissingNumberByTaxon(i)/gt.getSiteNumber());
            hetTaxon.add(gt.getHeterozygousProportionByTaxon(i));
        }


        String siteMissingPdf = new File (outDirS, genomeType+"_site_missing.pdf").getAbsolutePath();
        String siteHetPdf = new File (outDirS, genomeType+"_site_het.pdf").getAbsolutePath();
        String taxaMissingPdf = new File (outDirS, genomeType+"_taxa_missing.pdf").getAbsolutePath();
        String taxaHetPdf = new File (outDirS, genomeType+"_taxa_het.pdf").getAbsolutePath();
        String mafPdf = new File (outDirS, genomeType+"__site_maf.pdf").getAbsolutePath();
        String taxaHetFileS = new File (outDirS, genomeType+"_site_het.txt").getAbsolutePath();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(taxaHetFileS);
            bw.write("Taxa\tHeterozygousProportion");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < gt.getTaxaNumber(); i++) {
                sb.setLength(0);
                sb.append(gt.getTaxonName(i)).append("\t").append(hetTaxon.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        Histogram h = new Histogram(missingSite.toArray());
        h.setTitle(genomeType+"_missing"+"_bySite");
        h.setXLab("Missing rate");
        h.setYLab("Proportion");
        h.saveGraph(siteMissingPdf);
        h = new Histogram(hetSite.toArray());
        h.setTitle(genomeType+"_heterozygous"+"_bySite");
        h.setXLab("Heterozygous proportion by site");
        h.setYLab("Proportion");
        h.saveGraph(siteHetPdf);
        h = new Histogram(maf.toArray());
        h.setTitle(genomeType+"_MAF"+"_bySite");
        h.setXLab("MAF");
        h.setYLab("Proportion");
        h.saveGraph(mafPdf);
        h = new Histogram(missingTaxon.toArray());
        h.setTitle(genomeType+"_missing"+"_byTaxa");
        h.setXLab("Missing rate");
        h.setYLab("Proportion");
        h.saveGraph(taxaMissingPdf);
        h = new Histogram(hetTaxon.toArray());
        h.setTitle(genomeType+"_heterozygous"+"_byTaxa");
        h.setXLab("Heterozygous proportion by taxa");
        h.setYLab("Proportion");
        h.saveGraph(taxaHetPdf);
    }

    public void filterVCF () {
        String abInVCFDirS = "/Volumes/VMap2_Fei/vcf/001_fastcall/ab";
        String abdInVCFDirS = "/Volumes/VMap2_Fei/vcf/001_fastcall/abd";
        String dInVCFDirS = "/Volumes/VMap2_Fei/vcf/001_fastcall/d";

        String abOutVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/ab";
        String abdOutVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/abd";
        String dOutVCFDirS = "/Volumes/VMap2_Fei/vcf/002_byDepth/d";

        String reliableDirS = "/Volumes/VMap2_Fei/reliableSites/ABD_intersect";

        this.filter(abInVCFDirS, abOutVCFDirS, reliableDirS);
        this.filter(abdInVCFDirS, abdOutVCFDirS, reliableDirS);
        this.filter(dInVCFDirS, dOutVCFDirS, reliableDirS);

    }

    private void filter (String InVCFDirS, String outVCFDirS, String reliableDirS) {
        List<File> fList  = IOUtils.getFileListInDirEndsWith(InVCFDirS, ".gz");
        fList.parallelStream().forEach(f -> {
            String outVCF = new File (outVCFDirS, f.getName()).getAbsolutePath();
            String reliableFileS = new File (reliableDirS, f.getName().split("\\.")[0]+"_intersect_reliable.txt.gz").getAbsolutePath();
            try {
                BitSet bs = new BitSet();
                BufferedReader br  = IOUtils.getTextGzipReader(reliableFileS);
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (Integer.parseInt(temp) == 1) {
                        bs.set(cnt);
                    }
                    cnt++;
                }
                br.close();
                System.out.println(reliableFileS);
                br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outVCF);
                while ((temp = br.readLine()).startsWith("##")) {
                    bw.write(temp);
                    bw.newLine();
                }
                bw.write(temp);
                bw.newLine();
                cnt = 0;
                List<String> l = new ArrayList<>();
                int index = 0;
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp.substring(0, 50));
                    index = Integer.parseInt(l.get(1))-1;
                    if (bs.get(index)) {
                        bw.write(temp);
                        bw.newLine();
                    }
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(cnt);
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
}
