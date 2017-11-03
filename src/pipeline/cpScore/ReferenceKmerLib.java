package pipeline.cpScore;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.set.hash.HashIntSet;
import com.koloboke.collect.set.hash.HashIntSets;
import format.Fasta;
import utils.BaseEncoder;

import java.util.Arrays;

public class ReferenceKmerLib {
    int kmerLength = 32;


    public ReferenceKmerLib(int kmerLength, String inputGenomeFileS, String libFileS) {
        this.kmerLength = kmerLength;
        this.createKmerSet(inputGenomeFileS);
    }

    void createKmerSet (String inputGenomeFileS) {
        inputGenomeFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/maize_chr12.fa";
        Fasta f = new Fasta(inputGenomeFileS);
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
        System.out.println("Building kmer list from reference...");
        System.out.println("KmerLength = "+String.valueOf(kmerLength)+ " bp");
        int totalLength = (int)f.getTotalSeqLength();
        HashIntSet kmerSet = HashIntSets.newMutableSet(totalLength);
        long start = System.nanoTime();
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                boolean flag = false;
                for (int j = i; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                byte[] subArray = Arrays.copyOfRange(bArray, i, i + kmerLength);
                int kmerL = BaseEncoder.getIntSeqFromByteArray(subArray);
                if (!kmerSet.contains(kmerL)) kmerSet.add(kmerL);
                int pos = i+1;
                if (pos%10000000 == 0) {
                    System.out.println("Chromosome: "+f.getName(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos));
                }
            }
        }
    }

}
