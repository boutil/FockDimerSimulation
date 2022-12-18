package de.jtem.riemann.schottky;

import de.jtem.mfc.field.Complex;

public class SchottkyDimersQuad extends SchottkyDimers {
    public double alphaMinus;
    public double betaMinus;
    public double alphaPlus;
    public double betaPlus;
    public Complex[] angles;

    public SchottkyDimersQuad(SchottkyData data, double[] angles) {
        super(data, angles);
        alphaMinus = angles[0]; 
        betaMinus = angles[1];
        alphaPlus = angles[2];
        betaPlus = angles[3];
        this.angles = new Complex[] {new Complex(angles[0], 0), new Complex(angles[1], 0), 
            new Complex(angles[2], 0), new Complex(angles[3], 0)};

        Complex P0 = new Complex(0, 1);

        amoebaMap = new AmoebaMap(this, P0);
    }

    @Override
    public Complex amoebaMap(Complex P) throws Exception {
        // P needs to be in fundamental domain.
        if (!isInFundamentalDomain(P)) {
            throw new Exception("P needs to be in fundamental domain");
        }
        Complex r = new Complex();
        amoebaMap.quadGrid(r, P, angles, acc);
        return r;
    }


}
