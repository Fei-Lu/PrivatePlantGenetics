/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import pgl.utils.Dyad;
import pgl.utils.IOUtils;

/**
 *
 * @author feilu
 */
public class VMapDBUtils {
    
    public static Dyad<String, List<String>> getDBInfo (String dbFileS) {
        String header = null;
        List<String> recordList = new ArrayList();
        BufferedReader br = null;
        try {
            if (dbFileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(dbFileS);
            }
            else if (dbFileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(dbFileS);
            }
            else {
                return null;
            }
            String temp = null;
            header = br.readLine();
            while ((temp = br.readLine()) != null) {
                recordList.add(temp);
            }
            Dyad<String, List<String>> two = new Dyad(header, recordList);
            System.out.println("DB file read from " + dbFileS);
            return two;
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }
    
    public static void writeDB (String header, List<String> recordList, String dbFileS) {
        BufferedWriter bw = null;
        try {
            if (dbFileS.endsWith(".txt")) {
                bw = IOUtils.getTextWriter(dbFileS);
            }
            else if (dbFileS.endsWith(".gz")) {
                bw = IOUtils.getTextGzipWriter(dbFileS);
            }
            else {
                System.out.println("File suffix issue. quit.");
                System.exit(1);
            }
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < recordList.size(); i++) {
                bw.write(recordList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("DB output to " + dbFileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static int getColumnIndexInHeader (String header, String columnName) {
        String[] temp = header.split("\t");
        for (int i = 0; i < temp.length; i++) {
            if (temp[i].equals(columnName)) return i;
        }
        return -1;
    }

}
