package de.jtem.riemann.schottky;

import java.util.List;

import de.jtem.mfc.field.Complex;

public class SchottkyDimersQuad extends SchottkyDimers {

    public SchottkyDimersQuad(SchottkyData data, double[][] angles) {
        super(data, angles);

        Complex P0 = new Complex(0, 1);

        amoebaMap = new AmoebaMap(this, P0);
    }

    @Override
    public Complex amoebaMap(Complex P) throws Exception {
        // P needs to be in fundamental domain.
        // if (!isInFundamentalDomain(P)) {
        //     throw new Exception("P needs to be in fundamental domain");
        // }
        return amoebaMap.amoebaMapQuadGrid(P, acc);
    }

    public Complex aztecMap(Complex P) throws Exception {
        return amoebaMap.aztecMap(P, acc);
    }

    public Complex aztecArcticCurve(Complex P) throws Exception {
        return amoebaMap.aztecMap(P, acc);
        // Complex curvePoint = null;
        // if (Math.abs(P.im) < 0.001) {
        //     curvePoint = amoebaMap.aztecArcticCurveReal(P, getAngles(), acc);
        // } else {
        //     // find which circle this point belongs to and then use aztecArcticCurve
        //     for (int i = 0; i < numGenerators; i++) {
        //         if(Math.abs(getCenterOfCircle(i, false).minus(P).abs() - getRadius(i)) < 0.01) {
        //             curvePoint = amoebaMap.aztecArcticCurve(P, getAngles(), getCenterOfCircle(i, false), acc);
        //         }
        //     }
        // }
        // if (curvePoint == null) {
        //     throw new Exception();
        // }
        // return curvePoint;
    }


}
