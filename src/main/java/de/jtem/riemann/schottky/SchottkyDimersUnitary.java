package de.jtem.riemann.schottky;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import inUtil.ComplexUtil;

public class SchottkyDimersUnitary extends SchottkyDimers{

    // angles here interpreted as the arg of unitary complex angles.

    public Complex slopeTranslation;
    
    public SchottkyDimersUnitary(SchottkyData data, double[][] angles, double[] boundaryResidues) {
        super(data, angles);

        // TODO: P0 should be chosen such that it's in the fundamental domain.
        amoebaMap = new AmoebaMapHex(this, chooseP0(), boundaryResidues);

        checkMCurveProp();
    }

    public SchottkyDimersUnitary(SchottkyData data, double[][] angles) {
        this(data, angles, new double[]{1, -1, 1, -1, 1, -1});
    }

    private void checkMCurveProp() {
        double eps = 0.01;
        for (int i = 0; i < numGenerators; i++) {
            if (!getA(i).equals(getB(i).invert().conjugate(), eps)){
                System.out.println("not a unitary M-curve");
            }
        }
    }

    protected double[] getAnglesInOrder() {
        double[] orderedAngles = new double[numAngles];
        for (int i = 0; i < angles.length; i++) {
            for (int j = 0; j < angles[i].length; j++) {
                orderedAngles[i + j] = angles[i][j];
            }
        }
        return orderedAngles;
    }

    protected void getSlopeTranslation() {
        // computes a translation to be applied to results from amoebaMap.getSlope() 
        // to make them all land in the positive triangle.
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

    @Override
    public Complex[] parametrizeInnerRealOval(int generatorI, int numPoints) {
        // Map to complex conjugated case, get the two circles of same size. Parametrize one of them and map back to unitary case.
        Moebius gen = getGenerator(generatorI);
        Moebius circleToReal = new Moebius(new Complex(0, -1), new Complex(1, 0), new Complex(1, 0), new Complex(0, -1));
        Moebius inv = circleToReal.invert();
        Moebius realGen = circleToReal.times(gen).times(inv);
        // Moebius realGen = inv.times(gen).times(circleToReal);
        Complex[] ev = realGen.getEigenValues();
        Complex[] fp = realGen.getFixPoints();
        Schottky schottky = new Schottky(new double[] {fp[0].re, fp[0].im, fp[1].re, fp[1].im, ev[0].re, 0});
        Complex center = schottky.getCenterOfCircle(0);
        double radius = schottky.getRadius(0);
        if (center.im < 0) {
            center.assignConjugate();
        }
        Complex circleToRealOfA1 = circleToReal.applyTo(getA(0));
        Complex circleToRealOfA1Inv = circleToRealOfA1.invert().neg();
        Complex cirlceToRealOfA2 = circleToReal.applyTo(getA(1));
        Complex originalFixPoint = inv.applyTo(fp[1]);
        Complex[] points = getPointArrayOnCircle(center, radius, numPoints);
        return Arrays.stream(points).map((x) -> inv.applyTo(x)).toArray(Complex[]::new);
    }

    // public Complex[] parametrizeInnerRealOval(int generatorI, int numPoints) {
    //     return getPointArrayOnCircle(getCenterOfCircle(generatorI), getRadius(generatorI) * 0.5, numPoints);
    // }

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

    public Complex getPointWithSlope(Complex slope) {
        // build a lattice of points and compute their slopes. Find the point that has the needed slope.
        int angleNum = 500;
        int rNum = 100;
        Complex[][] slopes = new Complex[angleNum][rNum];
        Complex closestPoint = new Complex();
        double minimalSlopeDiff = Double.MAX_VALUE;
        // calculate slopes on grid
        for (int i = 0; i < angleNum; i++) {
            for (int j = 0; j < rNum; j++) {
                double theta = 2 * Math.PI / angleNum * i;
                double r = (double)(j+1) / (rNum + 1);
                Complex p = new Complex(Math.cos(theta), Math.sin(theta)).times(r);
                if(!isInFundamentalDomain(p)) {
                    continue;
                }
                slopes[i][j] = amoebaMap.getSlope(p, acc);
                Complex pointSlope = amoebaMap.getSlope(p, acc);
                // double torusSize = Math.PI;
                // Complex slopeDiff = ComplexUtil.compMod(pointSlope.minus(slope), torusSize);
                // double dist = Math.max(Math.min(slopeDiff.re, torusSize - slopeDiff.re), Math.min(slopeDiff.im, torusSize - slopeDiff.im));
                // dist = Math.max(dist, Math.PI - Math.abs(slopeDiff.re - slopeDiff.im));
                double dist = pointSlope.minus(slope).abs();
                if (dist < minimalSlopeDiff) {
                    closestPoint = p;
                    minimalSlopeDiff = dist;
                }
            }
        }
        System.out.println("winding point: " + closestPoint + ". slopeDiff: " + minimalSlopeDiff);
        return closestPoint;
    }

    public SchottkyDimersDoubleCoverUnitary getDoubleCover() {
        // Here we need to first find the point P0 of winding where slope is 0 center of the Newton triangle. 
        // Then find transformation h which fixes the unit circle and maps P0 to 0.
        // Then apply sqrt() to get winding point at 0. This should be half of the new double cover.
        Complex[][] anglesC = getAngles();
        Complex slopeP0 = amoebaMap.getSlope(chooseP0(), acc);
        Complex pointBetweenGammaAlpha = anglesC[0][0].plus(anglesC[anglesC.length - 1][anglesC[anglesC.length - 1].length - 1]);
        pointBetweenGammaAlpha.assignDivide(pointBetweenGammaAlpha.abs() / 0.99999);
        Complex pointBetweenBetaGamma = anglesC[1][anglesC[1].length - 1].plus(anglesC[2][0]);
        pointBetweenBetaGamma.assignDivide(pointBetweenBetaGamma.abs() / 0.99999);
        Complex slopeGammaAlpha = amoebaMap.getSlope(pointBetweenGammaAlpha, acc);
        Complex slopeBetaGamma = amoebaMap.getSlope(pointBetweenBetaGamma, acc);
        Complex correctSlope = slopeP0.plus(slopeGammaAlpha).plus(slopeBetaGamma).divide(3);
        Complex correctWindingPoint = getPointWithSlope(correctSlope);

        // h maps P to 0 and unit circle to unit circle.
        Moebius h = new Moebius(new Complex(1, 0), correctWindingPoint.neg(), correctWindingPoint.neg().conjugate(), new Complex(1,0));

        Complex zero = h.applyTo(correctWindingPoint);

        double[][] newAngles = new double[angles.length * 2][];
        double h0 = h.applyTo(anglesC[0][0]).arg();
        for (int i = 0; i < angles.length; i++) {
            newAngles[i] = new double[angles[i].length];
            newAngles[i + angles.length] = new double[angles[i].length];
            for (int j = 0; j < angles[i].length; j++) {
                newAngles[i][j] = doubleMod(h.applyTo(anglesC[i][j]).arg() / 2, 2 * Math.PI);
                newAngles[i + angles.length][j] = doubleMod(newAngles[i][j] + Math.PI, 2 * Math.PI);
            }
        }
        SchottkyData data = new SchottkyData(numGenerators * 2);
        for (int i = 0; i < numGenerators; i++) {
            // double theta = doubleMod(h.applyTo(getA(i)).arg(), Math.PI * 2) / 2;
            // h.times(getGenerator(i)).getA();
            Complex newA = h.applyTo(getA(i)).sqrt();
            data.setA(i, newA);
            data.setB(i, newA.invert().conjugate());
            // What should mu be?
            data.setMu(i, getMu(i));
            data.setA(i + numGenerators, newA.neg());
            data.setB(i + numGenerators, newA.invert().conjugate().neg());
            data.setMu(i + numGenerators, getMu(i));
        }
        return new SchottkyDimersDoubleCoverUnitary(data, newAngles, amoebaMap.boundaryResidues);
    }

    public SchottkyDimersDoubleCoverUnitary getSimpleDoubleCover() {
        // Gets a double cover. Does not move original angles. They should already be in a single half-plane.
        Complex[][] anglesC = getAngles();
        double[][] newAngles = new double[angles.length * 2][];
        for (int i = 0; i < angles.length; i++) {
            newAngles[i] = new double[angles[i].length];
            newAngles[i + angles.length] = new double[angles[i].length];
            for (int j = 0; j < angles[i].length; j++) {
                newAngles[i][j] = doubleMod(anglesC[i][j].arg(), 2 * Math.PI);
                newAngles[i + angles.length][j] = doubleMod(newAngles[i][j] + Math.PI, 2 * Math.PI);
            }
        }
        SchottkyData data = new SchottkyData(numGenerators * 2);
        for (int i = 0; i < numGenerators; i++) {
            // double theta = doubleMod(h.applyTo(getA(i)).arg(), Math.PI * 2) / 2;
            // h.times(getGenerator(i)).getA();
            Complex newA = getA(i);
            data.setA(i, newA);
            data.setB(i, newA.invert().conjugate());
            // What should mu be?
            data.setMu(i, getMu(i));
            data.setA(i + numGenerators, newA.neg());
            data.setB(i + numGenerators, newA.invert().conjugate().neg());
            data.setMu(i + numGenerators, getMu(i));
        }
        return new SchottkyDimersDoubleCoverUnitary(data, newAngles, amoebaMap.boundaryResidues);
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

    protected double doubleMod(double x, double y) {
        return x - Math.floor(x/y) * y;
    }
}
