package unclassified;

import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;

public class Unclassified01 {

    public Unclassified01 () {

    }

    public void baiyuanzhang () {
        String infileS = "/Users/feilu/Documents/document_work/lab_management/research/funding/2020/白院长/2002-2007.txt";
        String outfileS = "/Users/feilu/Documents/document_work/lab_management/research/funding/2020/白院长/2002-2007_efficiency.txt";
        RowTable<String> t = new RowTable<>(infileS);
        double[][] nitrogen = new double[3][];
        double[][] pesti = new double[3][];
        double[][] wheat = new double[3][];
        double[][] rice = new double[3][];
        double[][] maize = new double[3][];
        double[][] soy = new double[3][];
        int[] years = t.getColumnAsIntArray(0);
        String[] species = {"Wheat", "Rice", "Maize", "Soybean"};
        String[] countries = {"China", "US", "World"};
        String[] cate = {"Nitrogen", "Pesticide"};

        for (int i = 0; i < 3; i++) {
            nitrogen[i] = t.getColumnAsDoubleArray(i+4);
            pesti[i] = t.getColumnAsDoubleArray(i+1);
            wheat[i] = t.getColumnAsDoubleArray(i+7);
            rice[i] = t.getColumnAsDoubleArray(i+10);
            maize[i] = t.getColumnAsDoubleArray(i+13);
            soy[i] = t.getColumnAsDoubleArray(i+16);
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            sb.append("Years");
            for (int i = 0; i < species.length; i++) {
                for (int j = 0; j < countries.length; j++) {
                    for (int k = 0; k < cate.length; k++) {
                        sb.append("\t").append(species[i]).append("_").append(countries[j]).append("_").append(cate[k]).append(" (t/kg)");
                    }
                }
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < years.length; i++) {
                sb.setLength(0);
                sb.append(years[i]);
                for (int j = 0; j < species.length; j++) {
                    double[][] speciesValue = null;
                    if (j == 0) speciesValue = wheat;
                    else if (j == 1) speciesValue = rice;
                    else if (j == 2) speciesValue = maize;
                    else if (j == 3) speciesValue = soy;
                    for (int k = 0; k < countries.length; k++) {
                        for (int l = 0; l < cate.length; l++) {
                            double value = 0;
                            if (l == 0) {
                                value = speciesValue[k][i]/nitrogen[k][i];
                            }
                            else if (l == 1) {
                                value = speciesValue[k][i]/pesti[k][i];
                            }
                            sb.append("\t").append((float)value);
                        }
                    }
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
