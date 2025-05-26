package pgl.tool.dev;

import pgl.graph.r.BarPlot;
import pgl.graph.r.ScatterPlot;

import java.util.Random;

public class GraphicsDev {
    int size = 1000;
    double[] x = new double[size];
    double[] y = new double[size];
    double[] z = new double[size];
    double r = 0.8;

    public GraphicsDev() {
        this.generateData();
        this.barPlot();

    }

    public void barPlot () {
        double[] value = {1,5,8};
        String[] name = {"a", "b", "c"};
        ScatterPlot p = new ScatterPlot (x, y);
        p.showGraph();
//        BarPlot p = new BarPlot(value, name);
//        p.setPlottingCharacter(3);
//        System.out.println(p.getPlotStatement());
//        p.showGraph();
    }

    /**
     * Generates x,y,z with a pearson's r between x-y, y-z
     */
    public void generateData () {
        Random random = new Random();

        for (int i = 0; i < x.length; i++) {
            x[i] = random.nextGaussian();
            y[i] = random.nextGaussian();
            z[i] = random.nextGaussian();
        }

        // Cholesky to create the correlation
        for (int i = 0; i < size; i++) {
            y[i] = r * x[i] + Math.sqrt(1 - r*r) * y[i];
            z[i] = r * y[i] + Math.sqrt(1 - r*r) * z[i];
        }
    }


}
