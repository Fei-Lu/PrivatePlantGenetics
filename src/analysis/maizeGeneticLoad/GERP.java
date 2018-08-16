/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maizeGeneticLoad;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import format.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.CrossMapUtils;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author feilu
 */
class GERP {
    
    public GERP () {
        this.convertToV4GerpFile();
    }
    
    public void convertToV4GerpFile () {
        //this.mkAGPV3BED();
        //this.crossMapConvert();
        this.makeV3V4Map();
        //this.mkV4Gerp();
    }
    
    public void mkV4Gerp () {
        String inGerpDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/gerp/agpV3";
        String v3v4MapDirS = "/Users/feilu/Documents/database/maize/crossMap/V3V4Map";
        String outGerpDirS = "/Users/feilu/Documents/analysisL/production/maizeLoad/gerp/agpV4";
        File[] fs = new File(v3v4MapDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "map");
        List<File> fList = Arrays.asList(fs);
        HashIntIntMap[] v3v4Maps = new HashIntIntMap[fList.size()];
        fList.parallelStream().forEach(f -> {
            int index = Integer.valueOf(f.getName().replaceFirst(".map", "").replaceFirst("chr", ""))-1;
            v3v4Maps[index] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    l = PStringUtils.fastSplit(temp);
                    if (l.get(2).startsWith("N")) {
                        v3v4Maps[index].put(Integer.parseInt(l.get(1)), -1);
                    }
                    else {
                        v3v4Maps[index].put(Integer.parseInt(l.get(1)), Integer.parseInt(l.get(2)));
                    }
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        fs = new File(inGerpDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "full");
        fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int index = Integer.valueOf(f.getName().replaceFirst(".msa.in.rates.full", "").replaceFirst("roast.chrom.", ""))-1;
            HashIntIntMap cMap = v3v4Maps[index];
            String outfileS = new File (outGerpDirS, "chr"+PStringUtils.getNDigitNumber(3, index+1)+"_Gerp.txt").getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Chr\tPos\tTreeLength\tValue");
                bw.newLine();
                int cnt = 0;
                String temp = null;
                List<String> l = null;
                int v4Pos = 0;
                StringBuilder sb = null;
                int chr = index + 1;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    v4Pos = cMap.get(cnt);
                    if (v4Pos == -1) continue;
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder();
                    sb.append(chr).append("\t").append(v4Pos).append("\t").append(l.get(0)).append("\t").append(l.get(1));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void makeV3V4Map () {
        String inDirS = "/Users/feilu/Documents/database/maize/crossMap/temp";
        String inCMDirS = "/Users/feilu/Documents/database/maize/crossMap/temp2";
        String outDirS = "/Users/feilu/Documents/database/maize/crossMap/V3V4Map";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "bed");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName().replaceFirst(".bed", ".map")).getAbsolutePath();
            String unmapFileS = new File(inCMDirS, f.getName()+".unmap").getAbsolutePath();
            String mapFileS = new File(inCMDirS, f.getName()).getAbsolutePath();
            TIntArrayList unPosList = new TIntArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(unmapFileS);
                String temp = null;
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    unPosList.add(Integer.parseInt(l.get(1)));
                }
                br.close();
                int[] unPos = unPosList.toArray();
                Arrays.sort(unPos);
                br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedReader br2 = IOUtils.getTextReader(mapFileS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("Chr_V3\tPos_V3\tChr_V4\tPos_V4");
                bw.newLine();
                String temp1 = null;
                StringBuilder sb = new StringBuilder();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder(l.get(0));
                    int query = Integer.parseInt(l.get(1));
                    int index = Arrays.binarySearch(unPos, query);
                    int pos = query+1;
                    sb.append("\t").append(pos).append("\t");
                    if (index < 0)  {
                        while (cnt < pos) {
                            cnt++;
                            temp1 = br2.readLine();
                            l = PStringUtils.fastSplit(temp1);
                        }
                        sb.append(l.get(0)).append("\t").append(Integer.parseInt(l.get(1))+1);
                    
                    }
                    else {
                        sb.append("NA").append("\t").append("NA");
                        cnt++;
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                br2.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void crossMapConvert () {
        String intDirS = "/Users/feilu/Documents/database/maize/crossMap/temp";
        String outDirS = "/Users/feilu/Documents/database/maize/crossMap/temp2";
        new File (outDirS).mkdir();
        File[] fs = new File(intDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "bed");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName()).getAbsolutePath();
            CrossMapUtils cm = new CrossMapUtils(f.getAbsolutePath(), outfileS);
            cm.setMaizeV3ToV4();
            cm.convert();
        });
    }
    
    public void mkAGPV3BED () {
        String infileS = "/Users/feilu/Documents/database/maize/infoFile/ChrLenCentPosi_agpV3.txt";
        String outDirS = "/Users/feilu/Documents/database/maize/crossMap/temp";
        RowTable<String> t = new RowTable(infileS);
        List<Integer> chrList = new ArrayList();
        int[] chrLength = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            String s = t.getCell(i, 0);
            chrList.add(Integer.parseInt(s));
            chrLength[i] = Integer.parseInt(t.getCell(i, 1));
        }
        chrList.parallelStream().forEach(chr -> {
            String outfileS = new File (outDirS, "chr"+PStringUtils.getNDigitNumber(3, chr)+".bed").getAbsolutePath();
            try {
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < chrLength[chr-1]; i++) {
                    sb = new StringBuilder();
                    sb.append(chr).append("\t").append(i).append("\t").append(i);
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
}
