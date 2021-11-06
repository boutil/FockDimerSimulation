package de.jtem.riemann.schottky;

import de.jtem.blas.ComplexMatrix;
import de.jtem.mfc.field.Complex;

public class SchottkyDimers extends Schottky{

    // We want to assume uniformization of type U2 (i.e. B_n = A_n.conj(), mu_n real)

    // Choice of special points on real line. These are the train track angles on the hexagonal grid.
    private double alpha, beta, gamma;

    // Object that computes the amoebaMap having access to SchottkyDimers object.
    AmoebaMap amoebaMap;

    // angles for now are angles alpha, beta, gamma in R.
    public SchottkyDimers(SchottkyData schottkyData, double[] angles) {
        super(schottkyData);
        alpha = angles[0];
        beta = angles[1];
        gamma = angles[2];
        Complex P0 = new Complex(0, 1);

        amoebaMap = new AmoebaMap(this, P0);
    }


    public Complex amoebaMapHexGrid(Complex P) throws Exception {
        // P needs to be in fundamental domain.
        if (!isInFundamentalDomain(P)) {
            throw new Exception("P needs to be in fundamental domain");
        }
        Complex r = new Complex();
        amoebaMap.hexGrid(r, P, new Complex(alpha, 0), new Complex(beta, 0), new Complex(gamma, 0), acc);
        return r;
    }

    public static void main(String[] args) {
        Complex A = new Complex(-15, 15);
        Complex B = new Complex(15, 15);
        Complex muA = new Complex(0.00001, 0);
        Complex muB = new Complex(0.00001, 0);
        SchottkyData schottkyData = new SchottkyData( 2 );
        schottkyData.setA( 0, A );
        schottkyData.setB( 0, A.conjugate() );
        schottkyData.setMu( 0, muA );
        schottkyData.setA( 1, B );
        schottkyData.setB( 1, B.conjugate() );
        schottkyData.setMu( 1, muB );
        System.out.println(schottkyData.isClassical());
        Schottky schottky = new Schottky(schottkyData);
        System.out.println(schottky.isDifferentialSeriesEvaluable());
        System.out.println(schottky.isIntegralSeriesEvaluable());
        // ComplexMatrix perMat = schottky.getPeriodMatrix();
        // System.out.println(perMat);

        double[] angles = {-30, 5, 47};
        SchottkyDimers dimers = new SchottkyDimers(schottkyData, angles);

        int length = dimers.getCenters().length;
        for (int i = 0; i < length; i++) {
            System.out.println("" + dimers.getCenterOfCircle(i) + ", " + dimers.getRadius(i));
        }

        Complex P = new Complex(0, 0);
        try {
            Complex r = dimers.amoebaMapHexGrid(P);
            System.out.println(r);
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
}
