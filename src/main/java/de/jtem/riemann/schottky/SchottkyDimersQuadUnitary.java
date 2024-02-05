package de.jtem.riemann.schottky;

import java.util.LinkedList;
import java.util.List;

import de.jtem.mfc.field.Complex;

public class SchottkyDimersQuadUnitary extends SchottkyDimersQuad {
    
    // angles here interpreted as e^{i pi angle}
    public SchottkyDimersQuadUnitary(SchottkyData data, double[][] angles, double[] boundaryResidues) {
        super(data, angles);

        // TODO: P0 should be chosen such that it's in the fundamental domain.
        amoebaMap = new AmoebaMapQuad(this, chooseP0(), boundaryResidues);
    }

    public SchottkyDimersQuadUnitary(SchottkyData data, double[][] angles) {
        this(data, angles, new double[]{1, -1, 1, -1, 1, -1});
    }

    @Override
    public Complex[][] getAngles() {
        Complex[][] complexAngles = new Complex[angles.length][];
        for (int i = 0; i < angles.length; i++) {
            Complex[] iAngles = new Complex[angles[i].length];
            for (int j = 0; j < angles[i].length; j++) {
                iAngles[j] = new Complex(Math.cos(angles[i][j]), Math.sin(angles[i][j]));
            }
            complexAngles[i] = iAngles;   
        }
        return complexAngles;
    }

    @Override
    public Complex[] getPointArrayOnRealLine(double a, double b, int numPoints, double[] excludingInterval){
        List<Complex> res = new LinkedList<Complex>();
        double radius = 0.999;
        if (b < a) {
            b = b + 2 * Math.PI;
        }
        for (double angle = a; angle < b; angle += (b-a)/numPoints) {
            Complex point = new Complex(radius * Math.cos(angle), radius * Math.sin(angle));
            if(isInFundamentalDomain(point)) {
                res.add(point);
            }                
        }
        return res.toArray(Complex[]::new);
    }

    public Complex[] parametrizeInnerRealOval(int generatorI, int numPoints) {
        return getPointArrayOnCircle(getCenterOfCircle(generatorI), getRadius(generatorI), numPoints);
    }

    @Override
    public Complex chooseP0() {
        Complex[][] anglesC = getAngles();
        Complex lastAlpha = anglesC[0][anglesC[0].length - 1];
        Complex firstBeta = anglesC[1][0];
        Complex sum = lastAlpha.plus(firstBeta);
        double factor = -1;
        if(firstBeta.divide(lastAlpha).arg() > 0) {
            factor = 1;
        }
        return sum.times(factor).divide(sum.abs());
    }

    @Override
    void updateCircles() {
        // In the unitary case the two circles are not of the same size.
        // We compute the explicit centers and radius of the inner circle. Outer circle radius is wrong but for now we do not care.
        for (int i = 0; i < numGenerators; i++) {
            Complex A = getA(i);
            Complex B = getB(i);
            double a = A.abs();
            double b = B.abs();
            if (a > 1) {
                System.out.println("A needs to be inside the unit circle.");
            }
            Complex Mu = getMu(i);
            double mu = Mu.abs();
            Complex angle = A.divide(A.abs());
            Complex P = angle.neg().times(1-mu).divide(a * mu - b);
            double r = Math.sqrt((a - mu * b) / (a * mu - b) + (Math.pow((1 - mu) / (a * mu - b), 2)));
            center[i][0] = P;
            center[i][1] = P.divide((P.abs() - r) * (P.abs() + r));
            radius[i] = r;
        }
    }

}
