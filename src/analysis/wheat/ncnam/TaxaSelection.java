/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.ncnam;

import pgl.infra.table.RowTable;
import pgl.infra.tree.Newick;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
class TaxaSelection {
    
    public TaxaSelection () {
        //this.mkClusteringInput();//deprecated
        //this.renameTaxaInGenotypeFile();
        this.test();
        
        //this.selectTaxaFromTree();
    }
    
    public void test () {
        String mdsInfileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/358Lines/taxa_358.mds.txt";
        String infoFileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/source/344line_info.txt";
        String treeXmlFileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/358Lines/taxa_358.xml";
        String treeNwkFileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/358Lines/taxa_358.nwk.txt";
        String outputMDSFileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/selected.MDS.pdf";
        String outputTreeFileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/selected.Tree.xml";
        String selectedTaxaFileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/selected.taxa.txt";
        int numS = 4;
        double minHeight = 0;
        double maxHeight = 1110;
        RowTable<String> t = new RowTable<>(infoFileS);
        //14 height
        List<String> sList = new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String s = t.getCell(i, 14);
            if (s.equals("")) continue;
            double value = Double.valueOf(s);
            if (value < minHeight) continue;
            if (value > maxHeight) continue;
            sList.add(t.getCell(i, 9));
        }
        List<String> allTaxaList = t.getColumn(9);
        t = new RowTable(mdsInfileS);
        allTaxaList.retainAll(t.getColumn(0));
        sList.retainAll(t.getColumn(0));
        Collections.sort(sList);
        Collections.sort(allTaxaList);
        List<String> rList = new ArrayList();
        
        
        String nwkS = null;
        try {
            BufferedReader br = IOUtils.getTextReader(treeNwkFileS);
            nwkS = br.readLine();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        Newick nwk = new Newick (nwkS);
        for (int i = 0; i < allTaxaList.size(); i++) {
            int index = Collections.binarySearch(sList, allTaxaList.get(i));
            if (index < 0) {
                rList.add(allTaxaList.get(i));
            }
        }
        for (int i = 0; i < rList.size(); i++) {
            nwk.deleteTaxonNode(rList.get(i));
        }
        List<String> treeSList = nwk.selectTaxaWithMaxDiversity(numS);
        Collections.sort(treeSList);
        sList = treeSList;
        Collections.sort(sList);
        
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(sList, t.getCell(i, 0));
            if (index < 0) continue;
            xList.add(Double.parseDouble(t.getCell(i, 1)));
            yList.add(Double.parseDouble(t.getCell(i, 2)));
        }
        double[] x = xList.toArray();
        double[] y = yList.toArray();
        ScatterPlot sp = new ScatterPlot(x,y);
        sp.setXLab("PC 1");
        sp.setYLab("PC 2");
        sp.setTitle("MDS of selected taxa from 358 lines");
        sp.setXLim(-0.3, 0.3);
        sp.setYLim(-0.25, 0.15);
        sp.saveGraph(outputMDSFileS);
        String addS = "<property ref=\"style:font_color\" datatype=\"xsd:token\" applies_to=\"node\">#ff0000</property>";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(selectedTaxaFileS);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < sList.size(); i++) {
                bw.write(sList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            BufferedReader br = IOUtils.getTextReader(treeXmlFileS);
            bw = IOUtils.getTextWriter(outputTreeFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.contains("<name>")) {
                    String currentName = temp.split(">")[1];
                    currentName = currentName.split("<")[0];
                    bw.write(temp);
                    bw.newLine();
                    bw.write(br.readLine());
                    bw.newLine();
                    if (Collections.binarySearch(sList, currentName) < 0) {
                        
                    }
                    else {
                        bw.write(addS);
                        bw.newLine();
                    } 
                }
                else {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void selectTaxaFromTree () {
        String inputFileS = "/Users/feilu/Documents/analysisH/ncnam/taxaSelection/taxa_358.nwk.txt";
        String nwkS = null;
        int n = 10;
        try {
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            nwkS = br.readLine();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String temp = "((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154)";
            
        Newick nwk = new Newick (nwkS);
        List<String> taxaList = nwk.selectTaxaWithMaxDiversity(n);
        int a = 3;
    }
    
    public void renameTaxaInGenotypeFile () {
        String infileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/source/AA358.unique.mapped.ok.hmp.txt";
        String infoFileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/source/344line_info.txt";
        String outfileS = "/Users/feilu/Documents/analysisH/NC-NAM/taxaSelection/source/taxa_358.hmp.txt";
        RowTable<String> t = new RowTable(infoFileS);
        HashMap<String, String> taxaNameMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaNameMap.put(t.getCell(i, 10), t.getCell(i, 9));
        }
        t = new RowTable (infileS);
        List<String> header = t.getHeader();
        for (int i = 0; i < header.size(); i++) {
            String c = header.get(i);
            c = taxaNameMap.get(c);
            if (c == null) continue;
            header.set(i, c);
        }
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    public void mkClusteringInput () {
        String inputFileS = "/Users/feilu/Documents/analysisH/ncnam/taxaSelection/source/taxa_358.unique.mapped.ok.hmp.txt";
        String outputFileS = "/Users/feilu/Documents/analysisH/ncnam/taxaSelection/AA358.sub.wgcna.txt";
        int markerNum = 10000;
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            RowTable<String> t = new RowTable<>(inputFileS);
            List<String> header  = t.getHeader();
            List<String> taxaList = new ArrayList(header);
            for (int i = 0; i < 11; i++) {
                taxaList.remove(0);
            }
            StringBuilder sb = new StringBuilder("Marker");
            for (int i = 0; i < taxaList.size(); i++) {
                sb.append("\t").append(taxaList.get(i));
            }
            bw.write(sb.toString());
            bw.newLine();
            int nMarker = t.getRowNumber()-1;
            int step = nMarker/markerNum;
            
            for (int i = 0; i < markerNum; i++) {
                int cIndex = i*step+1;
                String[] temp = t.getCell(cIndex, 1).split("/");
                String homo1 = temp[0]+temp[0];
                String homo2 = temp[1]+temp[1];
                String het = temp[0]+temp[1];
                sb = new StringBuilder("m"+String.valueOf(i+1));
                for (int j = 0; j < taxaList.size(); j++) {
                    sb.append("\t");
                    String cGeno = t.getCell(cIndex, j+11);
                    if (cGeno.equals(homo1)) {
                        sb.append(0);
                    }
                    else if (cGeno.equals(homo2)) {
                        sb.append(2);
                    }
                    else if (cGeno.equals(het)){
                        sb.append(1);
                    }
                    else sb.append(1);
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
    
}
