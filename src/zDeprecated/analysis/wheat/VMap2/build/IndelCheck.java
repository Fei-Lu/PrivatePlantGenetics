package zDeprecated.analysis.wheat.VMap2.build;

import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;

import java.util.Arrays;


public class IndelCheck {

    public IndelCheck () {
        this.test();
    }

    public void test () {
        String hapInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/indel/chr036.vcf.gz";
        String fastCallInfileS = "/Users/feilu/Documents/analysisL/production/vmap2/indel/chr036.ABDgenome.vcf.gz";
        GenotypeGrid g1 = new GenotypeGrid(hapInfileS, GenoIOFormat.VCF_GZ);
        TIntArrayList positionList = new TIntArrayList();
        for (int i = 0; i < g1.getSiteNumber(); i++) {
            positionList.add(g1.getPosition(i));
        }
        int[] positions = positionList.toArray();

        GenotypeGrid g2 = new GenotypeGrid(fastCallInfileS, GenoIOFormat.VCF_GZ);
        int siteNumber = g2.getSiteNumber();
        System.out.println(siteNumber);
        int cnt = 0;
        int cntOccu = 0;
        for (int i = 0; i < g2.getSiteNumber(); i++) {
            int index = Arrays.binarySearch(positions, g2.getPosition(i));
            if (index < 0) continue;
            char alt = g2.getAlternativeAlleleBase(i);
            if (alt == 'D' || alt == 'I') {
                cnt++;
                if (g2.getAlternativeAlleleOccurrenceBySite(i) < 2) continue;
                cntOccu++;
            }
        }
        System.out.println(cnt);
        System.out.println(cntOccu);
    }
}
