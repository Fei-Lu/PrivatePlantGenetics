/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.data;

import pgl.format.dna.FastaBit;
import pgl.format.table.RowTable;
import java.io.File;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import pgl.utils.IOFileFormat;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

/**
 *
 * @author feilu
 */
class WheatReferenceGenome {
    
    public WheatReferenceGenome () {
        //this.checkIUPAC();
        //this.splitGenomeByChr();
        //this.splitGenome();
    }
    
    public void splitGenome () {
        String infileDirS = "/Users/feilu/Documents/database/wheat/reference/v1.0/byChr" ;
        String abdFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/ABD/abd_iwgscV1.fa.gz";
        String aFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/A/a_iwgscV1.fa.gz";
        String bFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/B/b_iwgscV1.fa.gz";
        String abFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/AB/ab_iwgscV1.fa.gz";
        String dFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/D/d_iwgscV1.fa.gz";
        String mapFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/chrNameMap.txt";
        RowTable<String> t = new RowTable<>(mapFileS);
        HashMap<String, String> currentOriMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            currentOriMap.put(t.getCell(i, 1), t.getCell(i, 0));
        }
        File[] fArray = new File(infileDirS).listFiles();
        fArray = IOUtils.listFilesEndsWith(fArray, ".gz");
        List<File> fList = Arrays.asList(fArray);
        FastaBit[] faArray = new FastaBit[fList.size()];
        AtomicInteger atoCnt = new AtomicInteger();
        fList.stream().forEach(f -> {
            try {
                FastaBit fb = new FastaBit(f.getAbsolutePath());
                faArray[atoCnt.getAndIncrement()] = fb;
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        FastaBit afb = new FastaBit (faArray);
        //afb.sortByName();
        afb.writeFasta(abdFileS, IOFileFormat.TextGzip);
        
        boolean[] ifOut = new boolean[afb.getSeqNumber()];
        
        for (int i = 0; i < afb.getSeqNumber(); i++) {
            String ori = currentOriMap.get(afb.getName(i));
            if (ori.split("_")[0].endsWith("A")) ifOut[i] = true;
            if (ori.contains("Mit")) ifOut[i] = true;
            if (ori.contains("Chl")) ifOut[i] = true;
            if (ori.contains("Un")) ifOut[i] = true;
        }
        afb.writeFasta(aFileS, ifOut, IOFileFormat.TextGzip);
        
        ifOut = new boolean[afb.getSeqNumber()];
        for (int i = 0; i < afb.getSeqNumber(); i++) {
            String ori = currentOriMap.get(afb.getName(i));
            if (ori.split("_")[0].endsWith("A") || ori.split("_")[0].endsWith("B")) ifOut[i] = true;
            if (ori.contains("Mit")) ifOut[i] = true;
            if (ori.contains("Chl")) ifOut[i] = true;
            if (ori.contains("Un")) ifOut[i] = true;
        }
        afb.writeFasta(abFileS, ifOut, IOFileFormat.TextGzip);

        ifOut = new boolean[afb.getSeqNumber()];
        for (int i = 0; i < afb.getSeqNumber(); i++) {
            String ori = currentOriMap.get(afb.getName(i));
            if (ori.split("_")[0].endsWith("D")) ifOut[i] = true;
            if (ori.contains("Mit")) ifOut[i] = true;
            if (ori.contains("Chl")) ifOut[i] = true;
            if (ori.contains("Un")) ifOut[i] = true;
        }
        afb.writeFasta(dFileS, ifOut, IOFileFormat.TextGzip);
        
        ifOut = new boolean[afb.getSeqNumber()];
        for (int i = 0; i < afb.getSeqNumber(); i++) {
            String ori = currentOriMap.get(afb.getName(i));
            if (ori.split("_")[0].endsWith("B")) ifOut[i] = true;
            if (ori.contains("Mit")) ifOut[i] = true;
            if (ori.contains("Chl")) ifOut[i] = true;
            if (ori.contains("Un")) ifOut[i] = true;
        }
        afb.writeFasta(bFileS, ifOut, IOFileFormat.TextGzip);
    }
    
    public void splitGenomeByChr () {
        String infileS = "/Users/feilu/Documents/database/wheat/reference/download/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta";
        String outDirS = "/Users/feilu/Documents/database/wheat/reference/v1.0/byChr"; 
        String mapFileS = "/Users/feilu/Documents/database/wheat/reference/v1.0/chrNameMap.txt";
        RowTable<String> t = new RowTable<>(mapFileS);
        HashMap<String, String> oriCurrentMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            oriCurrentMap.put(t.getCell(i, 0), t.getCell(i, 1));
        }
        FastaBit fn = new FastaBit(infileS);
        for (int i = 0; i < fn.getSeqNumber(); i++) {
            String name = fn.getName(i);
            name = oriCurrentMap.get(name);
            fn.setName(name, i);
            String outfileS = new File (outDirS, "chr"+PStringUtils.getNDigitNumber(3, Integer.valueOf(name))+".fa.gz").getAbsolutePath();
            fn.writeFasta(outfileS, i, IOFileFormat.TextGzip);
        }
    }
    
    public void checkIUPAC () {
        String infileS = "/Users/feilu/Documents/database/wheat/reference/download/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta";
        FastaBit f = new FastaBit(infileS);
        System.out.println(f.isThereN());
        String[] names = f.getNames();
        for (int i = 0; i < names.length; i++) {
            System.out.println(names[i]);
        }
    }
    
    
}
