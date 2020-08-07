package analysis.wheat.sv;

import org.ietf.jgss.GSSName;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeExport;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeTable;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.PStringUtils;


import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class VCFDeletion {

    public VCFDeletion () {
        this.test();
    }

    public void test () {
        //BitSet[][] bArray, GridDirection gd, String[] taxa, BiSNP[] snps)
        String infileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/genotype/source/chr001_exon_vmap2.1.vcf.gz";
        String binOutFileS = "/Users/feilu/Documents/analysisL/pipelineTest/svVCFDeletion/chr001_exon_vmap2.1_deletion.vcf.gz";
        GenotypeTable gt = new GenotypeGrid(infileS, GenoIOFormat.VCF_GZ);
        int siteNumber = gt.getSiteNumber();
        int taxaNumber = gt.getTaxaNumber();
        ArrayList<BitSet[]> dGenoList = new ArrayList<>();
        ArrayList<BiSNP> dSNPList = new ArrayList<>();
        int cnt = 0;
        for (int i = 0; i < siteNumber; i++) {
            if (gt.getMissingNumberBySite(i) == 0) continue;
            BiSNP snp = new BiSNP(gt.getChromosome(i), gt.getPosition(i), gt.getReferenceAlleleBase(i), 'D', null);
            dSNPList.add(snp);
            BitSet[] genoSite = new BitSet[3];
            genoSite[2] = new BitSet(taxaNumber);
            genoSite[1] = new BitSet(taxaNumber);
            genoSite[0] = new BitSet(taxaNumber);
            for (int j = 0; j < taxaNumber; j++) {
                if (gt.isMissing(i,j)) {
                    genoSite[0].set(j);
                }
            }
            genoSite[1] = (BitSet)genoSite[0].clone();
            dGenoList.add(genoSite);
            cnt++;
        }
        BitSet[][] bArray = dGenoList.toArray(new BitSet[dGenoList.size()][]);
        BiSNP[] snps = dSNPList.toArray(new BiSNP[dSNPList.size()]);
        String[] taxa = gt.getTaxaNames();
        GenotypeGrid gg = new GenotypeGrid(bArray, GenotypeGrid.GridDirection.BySite, taxa, snps);
        System.out.println(gg.getSiteNumber());
        GenotypeExport.output(gg, binOutFileS, GenoIOFormat.VCF_GZ);

    }
}
