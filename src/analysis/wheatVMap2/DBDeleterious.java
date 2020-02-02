/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatVMap2;

import pgl.format.genomeAnnotation.GeneFeature;
import pgl.format.range.Range;
import pgl.format.table.ColumnTable;
import pgl.format.table.RowTable;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.graphcis.tablesaw.TablesawUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Table;
import pgl.utils.Dyad;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;
/**
 *
 * @author feilu
 */
public class DBDeleterious {
    
    public DBDeleterious() {
//       this.extractInfoFromVMap2();
//       this.mkExonVCF();
//       this.mkExonAnnotation();
//       this.addSift();
//       this.addAncestral();
//       this.addDAF();
//       this.addGerp();
////    this.addPhyloP();
//      this.addRecombination();

    }

    public void addRecombination () {
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/test";
        String recombinationFileS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/005_recombination/iwgsc_refseqv1.0_recombination_rate_chrID.txt";
        ColumnTable<String> t = new ColumnTable<>(recombinationFileS);
        int chrNum = Integer.parseInt(t.getCell(t.getRowNumber()-1, 0));
        TIntArrayList[] startLists = new TIntArrayList[chrNum];
        TIntArrayList[] endLists = new TIntArrayList[chrNum];
        TFloatArrayList[] crossLists = new TFloatArrayList[chrNum];
        for (int i = 0; i < startLists.length; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            crossLists[i] = new TFloatArrayList();
        }
        int index = -1;
        for (int i = 0; i < t.getRowNumber(); i++) {
            index = Integer.parseInt(t.getCell(i, 0))-1;
            startLists[index].add(Integer.parseInt(t.getCell(i, 1)));
            endLists[index].add(Integer.parseInt(t.getCell(i, 2)));
            crossLists[index].add(Float.parseFloat(t.getCell(i, 3)));
        }
        List<File> fList = IOUtils.getFileListInDirEndsWith(dirS, ".gz");
        fList.parallelStream().forEach(f -> {
            //String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath());
            String header = two.getFirstElement();
            List<String> recordList = two.getSecondElement();
            String[] tem = header.split("\t");
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(f.getAbsolutePath());
                StringBuilder sb = new StringBuilder(header);
                sb.append("\t").append("RecombinationRate");
                bw.write(sb.toString());
                bw.newLine();
                int chrIndex = -1;
                int posIndex = -1;
                int currentPos = -1;
                List<String> l  = null;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(recordList.get(i));
                    chrIndex = Integer.parseInt(l.get(1))-1;
                    currentPos = Integer.parseInt(l.get(2));
                    posIndex = startLists[chrIndex].binarySearch(currentPos);
                    if (posIndex < 0) posIndex = -posIndex-2;
                    if (posIndex < 0) {
                        sb.append(recordList.get(i)).append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    if (currentPos < endLists[chrIndex].get(posIndex)) {
                        sb.append(recordList.get(i)).append("\t").append(crossLists[chrIndex].get(posIndex));
                    }
                    else {
                        sb.append(recordList.get(i)).append("\t").append("NA");
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
        });
    }

    /**
     * @deprecated
     */
    public void addPhyloP () {
        String phyloPDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/004_phylop/byChr";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/002_exonSNPAnnotation/";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String phyloPFileS = f.getName().split("_")[0]+"_phyloP.txt.gz";
            phyloPFileS = new File (phyloPDirS, phyloPFileS).getAbsolutePath();
            String header = null;
            List<String> recordList = new ArrayList();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                header = br.readLine();
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tPhyloP");
                header = sb.toString();
                String temp = null;
                TIntArrayList posList = new TIntArrayList();
                List<String> l = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                    l = PStringUtils.fastSplit(temp);
                    posList.add(Integer.parseInt(l.get(2)));
                }
                br.close();
                br = IOUtils.getTextGzipReader(phyloPFileS);
                br.readLine();
                int pos = -1;
                int index = -1;
                String[] phyloP = new String[posList.size()];
                for (int i = 0; i < phyloP.length; i++) phyloP[i] = "NA";
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    pos = Integer.parseInt(l.get(1));
                    index = posList.binarySearch(pos);
                    if (index < 0) continue;
                    phyloP[index] = l.get(2);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextGzipWriter(f.getAbsolutePath());
                bw.write(header);
                bw.newLine();
                for (int i = 0; i < posList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t").append(phyloP[i]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(f.getAbsolutePath());
        });
    }
    
    public void addGerp () {
        String gerpDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/003_gerp/byChr_26way";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation/";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String gerpFileS = f.getName().split("_")[0]+"_gerp.txt.gz";
            gerpFileS = new File (gerpDirS, gerpFileS).getAbsolutePath();
            String header = null;
            List<String> recordList = new ArrayList();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                header = br.readLine();
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tGerp");
                header = sb.toString();
                String temp = null;
                TIntArrayList posList = new TIntArrayList();
                List<String> l = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                    l = PStringUtils.fastSplit(temp);
                    posList.add(Integer.parseInt(l.get(2)));
                }
                br.close();
                br = IOUtils.getTextGzipReader(gerpFileS);
                br.readLine();
                int pos = -1;
                int index = -1;
                String[] gerp = new String[posList.size()];
                for (int i = 0; i < gerp.length; i++) gerp[i] = "NA";
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    pos = Integer.parseInt(l.get(1));
                    index = posList.binarySearch(pos);
                    if (index < 0) continue;
                    gerp[index] = l.get(2);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextGzipWriter(f.getAbsolutePath());
                bw.write(header);
                bw.newLine();
                for (int i = 0; i < posList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t").append(gerp[i]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(f.getAbsolutePath());
        });
    }
    
    public void addDAF () {
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation/";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String header = null;
            List<String> recordList = new ArrayList();
            String tem = null;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                header = br.readLine();
                List<String> l = PStringUtils.fastSplit(header);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tDerived_SIFT\tDAF\tDAF_ABD\t").append(l.get(9).replaceFirst("AAF", "DAF"));
                header = sb.toString();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextGzipWriter(f.getAbsolutePath());
                bw.write(header);
                bw.newLine();
                float daf = -1;
                float dafABD = -1;
                float dafOther = -1;
                String subMajor = null;
                String subMinor = null;
                String ancestral = null;
                String derivedSIFT = null;
                double subMaf = -1;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t");
                    l = PStringUtils.fastSplit(recordList.get(i));
                    ancestral = l.get(15);
                    if (ancestral.equals(l.get(3))) {
                        derivedSIFT = l.get(13);
                    }
                    else if (ancestral.equals(l.get(4))) {
                        derivedSIFT = l.get(14);
                    }
                    else derivedSIFT = "NA";
                    sb.append(derivedSIFT).append("\t");
                    if (l.get(5).equals(ancestral)) {
                        sb.append((float)Double.parseDouble(l.get(7))).append("\t");
                    }
                    else if (l.get(6).equals(ancestral)) {
                        sb.append((float)(1- Double.parseDouble(l.get(7)))).append("\t");
                    }
                    else sb.append("NA\t");
                    if (Double.parseDouble(l.get(8)) < 0.5) {
                        subMajor = l.get(3);
                        subMinor = l.get(4);
                        subMaf = Double.parseDouble(l.get(8));
                    }
                    else {
                        subMajor = l.get(4);
                        subMinor = l.get(3);
                        subMaf = 1 - Double.parseDouble(l.get(8));
                    }
                    if (ancestral.equals(subMajor)) {
                        sb.append((float)subMaf).append("\t");
                    }
                    else if (ancestral.equals(subMinor)) {
                        sb.append((float)(1-subMaf)).append("\t");
                    }
                    else sb.append("NA\t");
                    
                    tem = recordList.get(i);
                    if (Double.parseDouble(l.get(9)) < 0.5) {
                        subMajor = l.get(3);
                        subMinor = l.get(4);
                        subMaf = Double.parseDouble(l.get(9));
                    }
                    else {
                        subMajor = l.get(4);
                        subMinor = l.get(3);
                        subMaf = 1 - Double.parseDouble(l.get(9));
                    }
                    if (ancestral.equals(subMajor)) {
                        sb.append((float)subMaf);
                    }
                    else if (ancestral.equals(subMinor)) {
                        sb.append((float)(1-subMaf));
                    }
                    else sb.append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                System.out.println(tem);
                e.printStackTrace();
            }
        });
    }
    
    public void addAncestral () {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/002_ancestral/byChr";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation/";
        File[] fs = new File (inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String annoFileS = f.getName().split("_")[0]+"_SNP_anno.txt.gz";
            annoFileS = new File(dirS, annoFileS).getAbsolutePath();
            String header = null;
            List<String> recordList = new ArrayList();
            TIntArrayList posList = new TIntArrayList();
            String[] ancestral = null;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(annoFileS);
                header = br.readLine();
                String temp = null;
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                    l = PStringUtils.fastSplit(temp);
                    posList.add(Integer.parseInt(l.get(2)));
                }
                ancestral = new String[posList.size()];
                for (int i = 0; i < ancestral.length; i++) ancestral[i] = "NA";
                br.close();
                br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    int index = posList.binarySearch(pos);
                    if (index<0) continue;
                    ancestral[index] = l.get(2);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextGzipWriter(annoFileS);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tAncestral");
                bw.write(sb.toString());
                bw.newLine();
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t").append(ancestral[i]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    public void addSift() {
        String siftAltDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/001_sift/output_alt";
        String siftRefDirS = "/Users/feilu/Documents/analysisH/vmap2/003_annotation/001_sift/output_ref";
        String dirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        File[] fs = new File(siftAltDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".xls.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String dbFileS = f.getName().split("_")[0]+"_SNP_anno.txt.gz";
            dbFileS = new File (dirS, dbFileS).getAbsolutePath();
            String refFileS = new File (siftRefDirS, f.getName().split("_")[0]+"_exon_vmap2.1_reverseRefAlt_SIFTannotations.xls.gz").getAbsolutePath();
            RowTable<String> tAlt = new RowTable (f.getAbsolutePath());
            RowTable<String> tRef = new RowTable (refFileS);
            SIFTRecord[] records = new SIFTRecord[tAlt.getRowNumber()];
            for (int i = 0; i < records.length; i++) {
                SIFTRecord s = new SIFTRecord(Integer.parseInt(tAlt.getCell(i, 1)), tAlt.getCell(i, 3), tAlt.getCell(i, 4), tAlt.getCell(i,
                        7), tAlt.getCell(i, 8), tAlt.getCell(i, 12), tRef.getCell(i, 12));
                records[i] = s;
            }
            Arrays.sort(records);
            try {
                List<String> dbList = new ArrayList();
                String temp = null;
                BufferedReader br = IOUtils.getTextGzipReader(dbFileS);
                String header = br.readLine();
                while ((temp = br.readLine()) != null) {
                    dbList.add(temp);
                }
                br.close();
                BufferedWriter bw = IOUtils.getTextGzipWriter(dbFileS);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tRegion\tVariant_type\tAlt_SIFT\tRef_SIFT");
                bw.write(sb.toString());
                bw.newLine();
                List<String> l = null;
                for (int i = 0; i < dbList.size(); i++) {
                    l = PStringUtils.fastSplit(dbList.get(i));
                    SIFTRecord query = new SIFTRecord(Integer.parseInt(l.get(2)), l.get(4), l.get(10));
                    int index = Arrays.binarySearch(records, query);
                    if (index < 0) continue;
                    sb.setLength(0);
                    sb.append(dbList.get(i)).append("\t");
                    sb.append(records[index].region).append("\t");
                    sb.append(records[index].type).append("\t");
                    sb.append(records[index].altSift).append("\t").append(records[index].refSift);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
    
    class SIFTRecord implements Comparable<SIFTRecord> {

        public int pos;
        public String alt;
        public String transcript;
        public String region;
        public String type;
        public String altSift;
        public String refSift;

        
        public SIFTRecord(int pos, String alt, String transcript) {
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
        }
        
        public SIFTRecord(int pos, String alt, String transcript, String region, String type, String altSift, String refSift) {
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
            this.region = region;
            this.type = type;
            this.altSift = altSift;
            this.refSift = refSift;
        }

        @Override
        public int compareTo(SIFTRecord o) { //
            if (this.pos < o.pos) { 
                return -1;
                } else if (this.pos == o.pos) {
                    int index = this.alt.compareTo(o.alt);
                    if (index < 0) {
                        return -1;
                    }
                    else if (index > 0) {
                        return 1;
                    }
                    else return transcript.compareTo(o.transcript);
                } else {
                    return 1;
            }
        }
    }

    public void mkExonVCF () {
        String vmapDirS = "/Volumes/Fei_HDD_Mac/VMap2.1";
        String geneFeatureFileS = "/Users/feilu/Documents/database/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String hcGeneFileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        String outputDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/002_exonSNPVCF";
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();
        RowTable<String> t = new RowTable<>(hcGeneFileS);
        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(2));
        List<Integer> chrList = new ArrayList<>();
        for (int i = 0; i < chrSet.size(); i++) {
            chrList.add(i+1);
        }
        chrList.parallelStream().forEach(chrID -> {
            String inputVCF = new File (vmapDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_vmap2.1.vcf.gz").getAbsolutePath();
            String outputVCF = new File (outputDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_exon_vmap2.1.vcf.gz").getAbsolutePath();
            List<String> geneList = new ArrayList<>();
            List<String> tranList = new ArrayList<>();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int currentChr = Integer.parseInt(t.getCell(i, 2));
                if (currentChr < chrID) continue;
                else if (currentChr > chrID) break;
                geneList.add(t.getCell(i, 0));
                tranList.add(t.getCell(i, 1));
            }
            int geneIndex = -1;
            List<Range> allexonList = new ArrayList<>();
            for (int i = 0; i < geneList.size(); i++) {
                geneIndex = gf.getGeneIndex(geneList.get(i));
                for (int j = 0; j < gf.getTranscriptNumber(geneIndex); j++) {
                    if (!tranList.get(i).equals(gf.getTranscriptName(geneIndex, j))) continue;
                    List<Range> exonList = gf.getExonList(geneIndex, j);
                    allexonList.addAll(exonList);
                }
            }
            Collections.sort(allexonList);
            int[] starts = new int[allexonList.size()];
            int[] ends = new int[allexonList.size()];
            for (int i = 0; i < starts.length; i++) {
                starts[i] = allexonList.get(i).getRangeStart();
                ends[i] = allexonList.get(i).getRangeEnd();
            }
            try {
                BufferedReader br = IOUtils.getTextGzipReader(inputVCF);
                BufferedWriter bw = IOUtils.getTextGzipWriter(outputVCF);
                String temp = null;
                while ((temp = br.readLine()).startsWith("##")) {
                    bw.write(temp); bw.newLine();
                }
                bw.write(temp); bw.newLine();
                List<String> l = new ArrayList<>();
                int index = -1;
                int pos = -1;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp.substring(0, 60));
                    pos = Integer.parseInt(l.get(1));
                    index = Arrays.binarySearch(starts, pos);
                    if (index < 0) index = -index - 2;
                    if (index < 0) continue;
                    if (pos < ends[index]) {
                        bw.write(temp);bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(chrID+"  mkExonVCF");
        });
    }
    
    public void mkExonAnnotation() {
        String inDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/001_genicSNPByChr/";
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/003_exonSNPAnnotation";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = f.getName().split("_")[0]+"_SNP_anno.txt.gz";
            outfileS = new File (outDirS, outfileS).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
    
    public void extractInfoFromVMap2 () {
        int subLength = 150;
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/001_genicSNPByChr/";
        String vmapDirS = "/Volumes/Fei_HDD_Mac/VMap2.1/";
        File[] fs  = new File(vmapDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        String geneHCFileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        Table t = TablesawUtils.readTsv(geneHCFileS);
        System.out.println(t.structure());
        t.sortAscendingOn("Chr", "TranStart");
        IntColumn chrColumn = t.intColumn("chr");
        int chrNum = chrColumn.countUnique();
        TIntList[] startLists = new TIntList[chrNum];
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }
        for (int i = 0; i < t.rowCount(); i++) {
            startLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 3)));
            endLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 4)));
            tranLists[Integer.parseInt(t.getString(i, 2))-1].add(t.getString(i, 1));
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst(".vcf.gz", "_genicSNP.txt.gz")).getAbsolutePath();
            int[] dc = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(dc);
            StringBuilder sb = new StringBuilder();
            if (Arrays.binarySearch(dc, chrIndex+1) < 0) {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB\tTranscript");
            }
            else {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D\tTranscript");
            }
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()).startsWith("#")) {}
                
                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) {
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    posIndex = startLists[chrIndex].binarySearch(currentPos);
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    sb.append("\t").append(l.get(4)).append("\t");
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""),",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    }
                    else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]).append("\t").append(ll.get(7).split("=")[1]).append("\t").append(ll.get(8).split("=")[1]);
                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
}
