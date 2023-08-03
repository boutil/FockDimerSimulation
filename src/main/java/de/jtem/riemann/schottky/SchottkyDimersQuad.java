package de.jtem.riemann.schottky;

import java.util.List;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

import de.jtem.mfc.field.Complex;

public class SchottkyDimersQuad extends SchottkyDimers {

    public SchottkyDimersQuad(SchottkyData data, double[][] angles) {
        super(data, angles);

        Complex P0 = new Complex(0, 1);

        amoebaMap = new AmoebaMapQuad(this, P0);
    }

    public Complex aztecMap(Complex P) throws Exception {
        return amoebaMap.boundaryMap(P, acc);
    }

    public Complex aztecArcticCurve(Complex P) throws Exception {
         return amoebaMap.boundaryMap(P, acc);
    }

}
