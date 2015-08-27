/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author maubar
 */
public class MinowskiDistance {

    /**
     * Assumes x and y have the same length
     */
    public static Double minowskiDistance(double [] x, double [] y, double p){
        //Compute | xi - yi |
        double [] diffs = new double[x.length];
        for(int i = 0 ; i < x.length; i++ )
            diffs[i] = Math.abs(x[i] - y[i]);

        /**
         * Compute sum of powers
         * Math.pow seems safe about overflows,
         * returns Double.POSITIVE_INFINITY number is too big
         * (or neg inf if too small)
         */
        Double sum_of_pow = 0.0;
        for(Double d : diffs){
            sum_of_pow += Math.pow(d,p);
        }

        return Math.pow(sum_of_pow,(1/p));
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        double [] x = {1,0};
        double [] y = {0,1};
        System.out.println( minowskiDistance(x, y, 2));
    }

}
