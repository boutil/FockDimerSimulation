package de.jtem.riemann.schottky;


import de.jtem.mfc.field.Complex;

public class SchottkyDimersHex extends SchottkyDimers{

    // We want to assume uniformization of type U2 (i.e. B_n = A_n.conj(), mu_n real)

    // Choice of special points on real line. These are the train track angles on the hexagonal grid.
    private double alpha, beta, gamma;

    // Object that computes the amoebaMap having access to SchottkyDimers object.
    AmoebaMap amoebaMap;

    // angles for now are angles alpha, beta, gamma in R.
    public SchottkyDimersHex(SchottkyData schottkyData, double[][] angles) {
        super(schottkyData, angles);
        this.angles = angles;
    }

    @Override
    public Complex amoebaMap(Complex P) throws Exception {
        // P needs to be in fundamental domain.
        if (!isInFundamentalDomain(P)) {
            throw new Exception("P needs to be in fundamental domain");
        }
        Complex r = new Complex();
        amoebaMap.hexGrid(r, P, new Complex(alpha, 0), new Complex(beta, 0), new Complex(gamma, 0), acc);
        return r;
    }




}
