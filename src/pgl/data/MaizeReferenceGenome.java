package pgl.data;

import pgl.infra.dna.FastaBit;
import pgl.infra.dna.FastaByte;
import pgl.infra.dna.FastaRecordBit;
import pgl.infra.utils.IOFileFormat;

import java.util.ArrayList;

public class MaizeReferenceGenome {

    public MaizeReferenceGenome () {
//        this.checkIUPAC();
//        this.buildGenome();
//        this.printInfo();
    }

    public void printInfo () {
        String genomeFileS = "/Users/feilu/Documents/database/maize/reference/v5/genome/B73_v5_LuLab.fa.gz";
        FastaBit f = new FastaBit(genomeFileS);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            sb.append(f.getDescription(i)).append("\t").append(f.getSeqLength(i));
            System.out.println(sb.toString());
            sb.setLength(0);
        }
    }

    public void buildGenome () {
        String genomeInfileS = "/Users/feilu/Documents/database/maize/reference/v5/download/Zm-B73-REFERENCE-NAM-5.0.fa.gz";
        String mitoInfileS = "/Users/feilu/Documents/database/maize/reference/v5/download/mitochondrion.fasta.gz";
        String chloroInfileS = "/Users/feilu/Documents/database/maize/reference/v5/download/chloroplast.fasta.gz";
        String outfileS = "/Users/feilu/Documents/database/maize/reference/v5/genome/B73_v5_LuLab.fa.gz";
        ArrayList<FastaRecordBit> recordList = new ArrayList<FastaRecordBit>();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 60; i++) {
            sb.append("N");
        }
        String ns = sb.toString();
        sb.setLength(0);
        FastaBit f = new FastaBit(genomeInfileS);
        FastaRecordBit frb = null;
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String des = f.getDescription(i);
            if (des.startsWith("chr")) {
                des = des.replaceFirst("chr", "");
                frb = f.getFastaRecordBit(i);
                frb.setDescription(des);
                frb.setID(Integer.parseInt(des));
                recordList.add(frb);
            }
            else if (des.startsWith("scaf")) {
                if (i == f.getSeqNumber()-1) {
                    sb.append(f.getSeq(i));
                }
                else {
                    sb.append(f.getSeq(i)).append(ns);
                }
            }
            else {
                System.out.println("Error in the reference genome");
            }
        }
        System.out.println("Total Length: " + String.valueOf(f.getTotalSeqLength()));
        System.out.println("Total contigs " + String.valueOf(f.getSeqNumber()));
        frb = new FastaRecordBit("0", sb.toString(), 0);
        recordList.add(frb);
        f = new FastaBit(mitoInfileS);
        frb = f.getFastaRecordBit(0);
        frb.setDescription("11");
        frb.setID(11);
        recordList.add(frb);
        f = new FastaBit(chloroInfileS);
        frb = f.getFastaRecordBit(0);
        frb.setDescription("12");
        frb.setID(12);
        recordList.add(frb);
        f = new FastaBit(recordList.toArray(new FastaRecordBit[recordList.size()]));
        f.sortByDescriptionValue();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            System.out.println(f.getDescription(i) + " Length: " + String.valueOf(f.getSeqLength(i)));
        }
        f.writeFasta(outfileS, IOFileFormat.TextGzip);
    }

    public void checkIUPAC () {
        String infileS = "/Users/feilu/Documents/database/maize/reference/v5/download/Zm-B73-REFERENCE-NAM-5.0.fa.gz";
        FastaByte f = new FastaByte(infileS);
        System.out.println("Nuclear genome:");
        System.out.println("NonACGTN exists? : " + String.valueOf(f.isThereNonACGTNBase()));
        System.out.println("N exists? : " + String.valueOf(f.isThereN()));

        infileS = "/Users/feilu/Documents/database/maize/reference/v5/download/mitochondrion.fasta.gz";
        f = new FastaByte(infileS);
        System.out.println("NonACGTN exists? : " + String.valueOf(f.isThereNonACGTNBase()));
        System.out.println("N exists? : " + String.valueOf(f.isThereN()));

        infileS = "/Users/feilu/Documents/database/maize/reference/v5/download/chloroplast.fasta.gz";
        f = new FastaByte(infileS);
        System.out.println("Chloroplast: ");
        System.out.println("NonACGTN exists? : " + String.valueOf(f.isThereNonACGTNBase()));
        System.out.println("N exists? : " + String.valueOf(f.isThereN()));
    }
}
