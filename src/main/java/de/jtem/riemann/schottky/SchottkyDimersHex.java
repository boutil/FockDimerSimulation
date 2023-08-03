package de.jtem.riemann.schottky;


import de.jtem.mfc.field.Complex;

public class SchottkyDimersHex extends SchottkyDimers{

    // We want to assume uniformization of type U2 (i.e. B_n = A_n.conj(), mu_n real)

    // Object that computes the amoebaMap having access to SchottkyDimers object.
    AmoebaMap amoebaMap;

    // angles for now are angles alpha, beta, gamma in R.
    public SchottkyDimersHex(SchottkyData schottkyData, double[][] angles) {
        super(schottkyData, angles);
        Complex P0 = new Complex(0, 1);

        amoebaMap = new AmoebaMapHex(this, P0);
    }
}
