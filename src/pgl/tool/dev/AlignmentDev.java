package pgl.tool.dev;

import pgl.infra.align.g2.SAMPEAlignment;
import pgl.infra.align.g2.SAMSEAlignment;
import pgl.infra.align.g2.SEAlignRecord;

public class AlignmentDev {
    
    public AlignmentDev () {
        //this.bwaAlignment();
        //this.testSE();
        //this.testPE();
    }

    public void bwaAlignment () {
//        Single End
//        bwa mem -t 8 /Users/feilu/Documents/database/maize/reference/bwaLib/maizeAGPv4.fa /Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/rawFastq/K16BJS0001_1.fq.gz > /Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/sam/K16_se.sam
//        Paried End
//        bwa mem -t 8 /Users/feilu/Documents/database/maize/reference/bwaLib/maizeAGPv4.fa /Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/rawFastq/K16BJS0001_1.fq.gz /Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/rawFastq/K16BJS0001_2.fq.gz > /Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/sam/K16_pe.sam
    }

    public void testSE () {
        String singleFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/sam/K16_se.sam";
        SAMSEAlignment se = new SAMSEAlignment(singleFileS);
        int alignNumber = se.getAlignmentNumber();
        for (int i = 0; i < alignNumber; i++) {
            SEAlignRecord sar = se.getAlignmentRecord(i);
            System.out.println(sar.getAlignmentLength()+"\t"+sar.getAlignMatchNumber());
        }
    }
    
    public void testPE () {
        String peFileS = "/Users/feilu/Documents/analysisL/softwareTest/pgl/alignment/sam/K16_pe.sam";
        SAMPEAlignment pe = new SAMPEAlignment(peFileS);
        for (int i = 0; i < pe.getAlignmentNumber(); i++) {
            System.out.println(pe.getAlignmentRecord(i).getInsertSize());
        }
    }
}
