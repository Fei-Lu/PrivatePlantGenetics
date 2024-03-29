/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.genomeSize;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import com.koloboke.collect.map.hash.HashByteByteMap;
import pgl.infra.dna.BaseEncoder;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author feilu
 */
public class KmerCount {
    String taxon;
    byte kmerLength;
    long[] kmers;
    int[] counts;
    
    public KmerCount (String inputFileS, String taxon) {
        this.readFromText(inputFileS, taxon);
    }
    
    public void readFromText (String inputFileS, String taxon) {
        this.taxon = taxon;
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIBaseCodingMap();
        
        try {
            TLongArrayList kmerList = new TLongArrayList();
            TIntArrayList countList = new TIntArrayList();
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            String temp = null;
            List<String> temList = null;
            byte[] kmerByte = null;
            while ((temp = br.readLine()) != null) {
                temList = PStringUtils.fastSplit(temp);
                kmerLength = (byte)temList.get(0).length();
                kmerByte = temList.get(0).getBytes();
                for (int i = 0; i < kmerByte.length; i++) {
                    kmerByte[i] = ascIIByteMap.get(kmerByte[i]);
                }
                kmerList.add(BaseEncoder.getLongSeqFromBaseCodingArray(kmerByte));
                countList.add(Integer.parseInt(temList.get(1)));
                int a = 3;
            }
            br.close();
            kmers = kmerList.toArray();
            counts = countList.toArray();
            GenericSorting.quickSort(0, this.getKmerNumber(), compByKmer, swapper);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeText (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String s;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.getKmerNumber(); i++) {
                sb.setLength(0);
                s = BaseEncoder.getSequenceFromLong(this.getKmer(i)).substring(0, this.getKmerLength());
                sb.append(s).append("\t").append(this.getKmerCount(i));
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
    
    public void readBinaryFile (String infileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryReader(infileS);
            this.taxon = dis.readUTF();
            this.kmerLength = dis.readByte();
            int kmerNumber = dis.readInt();
            this.kmers = new long[kmerNumber];
            this.counts = new int[kmerNumber];
            for (int i = 0; i < kmerNumber; i++) {
                kmers[i] = dis.readLong();
                counts[i] = dis.readInt();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeBinaryFile (String outfileS) {
        try{
            DataOutputStream dos = IOUtils.getBinaryWriter(outfileS);
            dos.writeUTF(taxon);
            dos.writeByte(kmerLength);
            dos.write(this.getKmerNumber());
            for (int i = 0; i < this.getKmerNumber(); i++) {
                dos.writeLong(this.getKmer(i));
                dos.writeInt(this.getKmerCount(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public int getKmerIndex (long query) {
        return Arrays.binarySearch(this.kmers, query);
    }
    
    public String getTaxon () {
        return this.taxon;
    }
    
    public int getKmerNumber() {
        return this.kmers.length;
    }
    
    public int getKmerLength () {
        return this.kmerLength;
    }
    
    public long getKmer (int index) {
        return this.kmers[index];
    }
    
    public int getKmerCount (int index) {
        return this.counts[index];
        
    }
    
    protected Swapper swapper = new Swapper() {
        @Override
        public void swap(int a, int b) {
            long k = kmers[a];
            kmers[a] = kmers[b];
            kmers[b] = k;
            int c = counts[a];
            counts[a] = counts[b];
            counts[b] = c;
        }
    };
    
    protected IntComparator compByKmer = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            if (kmers[a] < kmers[b]) return -1;
            else if (kmers[a] > kmers[b]) return 1;
            else return 0;
        }
    };
    
    protected IntComparator compByCount = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            if (counts[a] < counts[b]) return -1;
            else if (counts[a] > counts[b]) return 1;
            else return 0;
        }
    };
    
}

