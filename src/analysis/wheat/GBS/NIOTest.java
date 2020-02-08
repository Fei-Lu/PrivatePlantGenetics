/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheat.GBS;

import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.Dyad;

/**
 *
 * @author feilu
 */
class NIOTest {
    
    public NIOTest () {
//        this.basicChannel();
//        this.readAndWrite();
        this.testBinaryWriteAndRead();
    }
    
    public void testBinaryWriteAndRead () {
        byte a = -1;
        int b = 3;
        String c = "abc";
        double d = 0.138d;
        String outfileS = "/Users/feilu/Documents/analysisL/production/nio/nio-data-b.txt";
        try {
            Dyad<FileChannel, ByteBuffer> iot = IOUtils.getNIOChannelBufferWriter(outfileS, 64);
            FileChannel fc = iot.getFirstElement();
            ByteBuffer bb = iot.getSecondElement();
            System.out.print(bb.position()+"\t");
            System.out.print(bb.limit()+"\t");
            System.out.println(bb.remaining()+"\n");
            
            bb.putInt(b);
            //bb.put(a);
            System.out.print(bb.position()+"\t");
            System.out.print(bb.limit()+"\t");
            System.out.println(bb.remaining()+"\n");
            
            bb.flip();
            bb.limit(bb.capacity());
            System.out.print(bb.position()+"\t");
            System.out.print(bb.limit()+"\t");
            System.out.println(bb.remaining()+"\n");
            
            fc.write(bb);
            System.out.print(bb.position()+"\t");
            System.out.print(bb.limit()+"\t");
            System.out.println(bb.remaining()+"\n");
            
            bb.clear();
            fc.close();
            
            Dyad<FileChannel, ByteBuffer> ioti = IOUtils.getNIOChannelBufferReader(outfileS, 128);
            FileChannel fci = ioti.getFirstElement();
            bb = ioti.getSecondElement();
            fci.read(bb);
            bb.flip();
            System.out.println(bb.getInt());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void readAndWrite () {
        try {
            RandomAccessFile inFile = new RandomAccessFile("/Users/feilu/Documents/analysisL/production/nio/nio-data.txt", "rw");
            FileChannel inChannel = inFile.getChannel();
            RandomAccessFile outFile = new RandomAccessFile("/Users/feilu/Documents/analysisL/production/nio/nio-data-out.txt", "rw");
            FileChannel outChannel = outFile.getChannel();
//            ByteBuffer buf = ByteBuffer.allocate(48);
            ByteBuffer buf = ByteBuffer.allocateDirect(48);
            //MappedByteBuffer buf = inChannel.map(FileChannel.MapMode.READ_WRITE, 0, 1024);
            int bytesRead = -1;
            while ((bytesRead = inChannel.read(buf)) != -1) {
                buf.flip();
                outChannel.write(buf);
                buf.clear();
            }
            inFile.close();
            outFile.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void basicChannel () {
        try {
            RandomAccessFile aFile = new RandomAccessFile("/Users/feilu/Documents/analysisL/production/nio/nio-data.txt", "rw");
            FileChannel inChannel = aFile.getChannel();

            ByteBuffer buf = ByteBuffer.allocate(48);

            int bytesRead = -1;
            while ((bytesRead = inChannel.read(buf)) != -1) {

                System.out.println("Read " + bytesRead);
                System.out.println(buf.limit());
                System.out.println(buf.position()+"\n");

                buf.flip();
                System.out.println(buf.limit());
                System.out.println(buf.position()+"\n");
                
                while(buf.hasRemaining()){
                    System.out.println((char) buf.get());
                    System.out.println(buf.limit());
                    System.out.println(buf.position()+"\n");
                }

                buf.clear();
                System.out.println(buf.limit());
                System.out.println(buf.position());

            }
            aFile.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
