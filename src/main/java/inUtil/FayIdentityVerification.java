package inUtil;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.blas.IntegerVector;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.Schottky;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.theta.Theta;
import de.jtem.riemann.theta.ThetaWithChar;

public class FayIdentityVerification {
    


    public static void main(String[] args) {

        // Create a Schottky uniformization with some chosen parameters.

        // Some Examples: uncomment only one of them:
        // M-curves with inversion in unit sphere as anti-holomorphic involution
        // -------------------------
        // // G1:
        // int genus = 1;
        // double theta1 = Math.PI;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.35);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParams = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000003, 0};
        // -------------------------
        // G2:
        int genus = 2;
        double theta1 = Math.PI;
        Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.35);
        Complex A1Refl = A1.invert().conjugate();
        double theta2 = 3 * Math.PI / 2;
        Complex A2 = new Complex(Math.cos(theta2), Math.sin(theta2)).times(0.35);
        Complex A2Refl = A2.invert().conjugate();
        double[] schottkyParams = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.0003, 0, A2.re, A2.im, A2Refl.re, A2Refl.im, 0.0003, 0};
        
        // -------------------------
        // // G1
        // Complex conjugation is anti-holomorphic involution
        // int genus = 1;
        // Complex A1 = new Complex(0,1);
        // Complex A1Conj = A1.conjugate();
        // double[] schottkyParams = {A1.re, A1.im, A1Conj.re, A1Conj.im, 0.01, 0};

        // --------------------------
        // Choose some 4 points to verify the Fay identity at.        
        // Complex[] points = new Complex[]{new Complex(-1, 0.5), new Complex(0, 0.3), new Complex(1, 0.5), new Complex(2, 0.5)};
        Complex[] points = new Complex[]{new Complex(-1, 0), new Complex(0, 0), new Complex(1, 0), new Complex(2, 0)};
        // Complex[] points = new Complex[]{new Complex(0.5, 0), new Complex(0.3, 0), new Complex(0, 0.5), new Complex(0.3, 0.3)};
        
        // Parameter choice over - rest is calculation
        // ====================================================


        Schottky schottky = new Schottky(new SchottkyData(schottkyParams));

        // Create a theta function object
        ComplexMatrix B = schottky.getPeriodMatrix();
        Theta theta = new Theta(schottky.getPeriodMatrix(), 1e-7, false);
        // Create a theta function with odd characteristic (1, 0, 0, ...).
        ThetaWithChar thetaWithChar = new ThetaWithChar(schottky.getPeriodMatrix(), 1e-7, false);
        int[] oddAlpha = new int[schottky.getNumGenerators()];
        oddAlpha[0] = 1;
        IntegerVector e1 = new IntegerVector(oddAlpha);
        thetaWithChar.setAlpha(e1);
        thetaWithChar.setBeta(e1);



        // fixed divisor
        ComplexVector D = new ComplexVector(genus);

        // Compute all Abel maps
        ComplexVector[] abelMaps = new ComplexVector[4];
        for (int i = 0; i < abelMaps.length; i++) {
            abelMaps[i] = new ComplexVector(genus);
            schottky.abelMap(abelMaps[i], points[i]);
        }

        // --------------------------
        // Compute the three terms and print the sum.
        Complex firstSummand = theta.theta(abelMaps[2].plus(abelMaps[3]).plus(D)).times(theta.theta(abelMaps[0].plus(abelMaps[1]).plus(D)));
        firstSummand.assignTimes(thetaWithChar.theta(abelMaps[3].minus(abelMaps[2])).times(thetaWithChar.theta(abelMaps[1].minus(abelMaps[0]))).divide(thetaWithChar.theta(abelMaps[2].minus(abelMaps[1]))).divide(thetaWithChar.theta(abelMaps[3].minus(abelMaps[0]))));

        Complex secondSummand = theta.theta(abelMaps[1].plus(abelMaps[3]).plus(D)).times(theta.theta(abelMaps[0].plus(abelMaps[2]).plus(D)));
        secondSummand.assignTimes(thetaWithChar.theta(abelMaps[1].minus(abelMaps[3])).times(thetaWithChar.theta(abelMaps[2].minus(abelMaps[0]))).divide(thetaWithChar.theta(abelMaps[2].minus(abelMaps[1]))).divide(thetaWithChar.theta(abelMaps[3].minus(abelMaps[0]))));

        Complex thirdSummand = theta.theta(abelMaps[1].plus(abelMaps[2]).plus(D)).times(theta.theta(abelMaps[0].plus(abelMaps[3]).plus(D)));

        Complex fay = firstSummand.plus(secondSummand).plus(thirdSummand);

        // Print the Fay sum. Hopefully it is close to zero.
        System.out.println("The sum of the Fay identity is:");
        System.out.println(fay);
    }

}
