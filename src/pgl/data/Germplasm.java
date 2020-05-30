package pgl.data;

import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;

class Germplasm {

    public Germplasm () {
        this.v1();
    }

    public void v1 () {
        String infileS = "/Users/feilu/Documents/document_work/lab_management/germplasm/wheat/source/germ_editByFei.txt";
        String outfileS = "/Users/feilu/Documents/document_work/lab_management/germplasm/wheat/germplasm_PGL.v1.txt";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while((temp = br.readLine()) != null) {
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
    }
}
