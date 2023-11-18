/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zDeprecated.analysis.wheat.VMap1;

/**
 *
 * @author feilu
 */
public class ConvergenceHyper {

    public ConvergenceHyper () {
        this.hyperTest();
    }

    public void hyperTest () {
        double as = 312;
        double ax = 3910;
        double ay = 2283;
        double g0 = 35345;
        System.out.println(this.getCHyper(as,ax,ay,g0));
    }

    public double getCHyper (double as, double ax, double ay, double g0) {
        double c = as-(ax*ay)/g0;
        double d = Math.sqrt(ax*ay*(g0-ax)*(g0-ay)/(g0*g0*(g0-1)));
        return c/d;
    }
}
