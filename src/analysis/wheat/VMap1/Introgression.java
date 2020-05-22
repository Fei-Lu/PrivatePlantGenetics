package analysis.wheat.VMap1;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.graph.r.DensityPlot;
import pgl.graph.r.Histogram;
import pgl.infra.range.Range;
import pgl.infra.range.Ranges;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

class Introgression {

    public Introgression () {
//        this.intervalSize();
        this.findIndividualWithMaxFd();
    }

    public void findIndividualWithMaxFd () {
        String inDirS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw";
        List<File> dirListA = IOUtils.getDirListInDirStartsWith(inDirS, "A");
        List<File> dirListB = IOUtils.getDirListInDirStartsWith(inDirS, "A");
        String sampleAFileS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/A025/A1.csv.gz";
        String sampleBFileS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/B001/AB1.csv.gz";
        Ranges ranA = this.getRanges(sampleAFileS, "A");
        Ranges ranB = this.getRanges(sampleBFileS, "B");
        String[] taxaA = this.getTaxaNames(dirListA.get(0));
        String[] taxaB = this.getTaxaNames(dirListA.get(0));

    }

    private String[] getTaxaNames (File dir) {
        List<File> fList = IOUtils.getFileListInDirEndsWith(dir.getAbsolutePath(), ".gz");
        String[] taxa = new String[fList.size()];
        for (int i = 0; i < taxa.length; i++) {
            taxa[i] = fList.get(i).getName().split("\\.")[0];
        }
        return taxa;
    }

    private Ranges getRanges (String sampleFileS, String subgenome) {
        RowTable<String> t = new RowTable<>(sampleFileS, ",");
        List<Range> rList = new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chr = Integer.parseInt(t.getCell(i,0));
            String type = RefV1Utils.getSubgenomeFromChrID(chr);
            if (type.equals(subgenome)) {
                Range r = new Range(chr, Integer.parseInt(t.getCell(i,1)), Integer.parseInt(t.getCell(i,2)));
                rList.add(r);
            }
        }
        return new Ranges(rList);
    }

    public void intervalSize() {
//        String file1 = "/Users/feilu/Desktop/untitled folder/A1.csv";
//        String file2 = "/Users/feilu/Desktop/untitled folder/A2.csv";
//        String file3 = "/Users/feilu/Desktop/untitled folder/AB1.csv";
//        String file4 = "/Users/feilu/Desktop/untitled folder/AB2.csv";
        String infileS = "/Volumes/Fei_HDD_Mac/VMap1.0/fd/all_individual/raw/B001/AB1.csv.gz";
        String intervalA = "/Users/feilu/Documents/analysisH/vmap1/fd/interval/A_interval_density.pdf";
        String intervalB = "/Users/feilu/Documents/analysisH/vmap1/fd/interval/B_interval_density.pdf";
        RowTable<String> t = new RowTable<>(infileS, ",");
        int[] chrs = t.getColumnAsIntArray(0);
        TIntHashSet cSet = new TIntHashSet(chrs);
        chrs = cSet.toArray();
        Arrays.sort(chrs);
        TDoubleArrayList AList = new TDoubleArrayList();
        TDoubleArrayList BList = new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String genomeType = RefV1Utils.getSubgenomeFromChrID(Integer.parseInt(t.getCell(i,0)));
            int dis = Integer.parseInt(t.getCell(i,2))-Integer.parseInt(t.getCell(i,1));
            if (genomeType.equals("A")) {
                AList.add(dis);
            }
            else if (genomeType.equals("B")) {
                BList.add(dis);
            }
        }
        DensityPlot d = new DensityPlot(AList.toArray());
        d.setSmoothN(5000);
        d.setTitle("Interval distribution of fd test in A subgenome");
        d.setXLab("Interval size (bp)");
        d.setYLab("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalA);
        d = new DensityPlot(BList.toArray());
        d.setSmoothN(5000);
        d.setTitle("Interval distribution of fd test in B subgenome");
        d.setXLab("Interval size (bp)");
        d.setYLab("Density");
        d.setXLim(0, 10000000);
        d.saveGraph(intervalB);

    }

}
