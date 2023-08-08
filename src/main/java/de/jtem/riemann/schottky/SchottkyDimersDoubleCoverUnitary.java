package de.jtem.riemann.schottky;

import java.util.LinkedList;
import java.util.List;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;

public class SchottkyDimersDoubleCoverUnitary extends SchottkyDimersUnitary{

    // angles here interpreted as the arg of unitary complex angles.

    private ComplexVector abelMapAtZero;
    
    public SchottkyDimersDoubleCoverUnitary(SchottkyData data, double[][] angles) {
        super(data, angles);

        Complex P0 = chooseP0();

        amoebaMap = new AmoebaMapHex(this, P0);

        ComplexVector abelMap0 = new ComplexVector(2, 0, 0);
        abelMap(abelMap0, new Complex(1,0));
        abelMapAtZero = new ComplexVector(1, abelMap0.getRe(0) + abelMap0.getRe(1),  abelMap0.getIm(0) + abelMap0.getIm(1));

        // It is assumed that the points in data and angles are doubled in z -> -z.
        checkMCurveProp();
        checkSlopeAtWinding();
    }

    private void checkMCurveProp() {
        double eps = 0.01;
        for (int i = 0; i < numGenerators; i++) {
            if (!getA(i).equals(getB(i).invert().conjugate(), eps)){
                System.out.println("not a unitary M-curve");
            }
        }
    }

    public void checkSlopeAtWinding() {
        Complex[][] anglesC = getAngles();
        Complex slope0 = amoebaMap.getSlope(new Complex(), acc);
        Complex slopeP0 = amoebaMap.getSlope(chooseP0(), acc);
        Complex pointBetweenGammaAlpha = anglesC[0][0].plus(anglesC[anglesC.length - 1][anglesC[anglesC.length - 1].length - 1]);
        pointBetweenGammaAlpha.assignDivide(pointBetweenGammaAlpha.abs());
        Complex pointBetweenBetaGamma = anglesC[1][anglesC[1].length - 1].plus(anglesC[2][0]);
        pointBetweenBetaGamma.assignDivide(pointBetweenBetaGamma.abs());
        System.out.println("slope at P0 is " + slopeP0);
        System.out.println("slope at betaGamma is " + amoebaMap.getSlope(pointBetweenBetaGamma, acc));
        System.out.println("slope at gammaAlpha is " + amoebaMap.getSlope(pointBetweenGammaAlpha, acc));
        System.out.println("slope at 0 should be: " + slopeP0.plus(amoebaMap.getSlope(pointBetweenGammaAlpha, acc)).plus(amoebaMap.getSlope(pointBetweenBetaGamma, acc)).divide(3));
        System.out.println("slope at 0 is " + slope0);
    }
    
    // For now we just assume genus 2. i.e. factor has genus 1.
    public ComplexMatrix getBMatrixOfFactor() {
        ComplexMatrix B = getPeriodMatrix();
        // Complex bPer = B.get(0, 0).plus(B.get(1, 0)).plus(B.get(0,1)).plus(B.get(1,1));
        Complex bPer = B.get(0, 0).plus(B.get(1, 0)).plus(B.get(0,1)).plus(B.get(1,1)).divide(2);
        return new ComplexMatrix(1, bPer);
    }

    public void abelMapOfFactor(ComplexVector v, Complex z) {
        ComplexVector map = new ComplexVector(2, 0, 0);
        abelMap(map, z);
        v.assign(map.getRe(0) + map.getRe(1), map.getIm(0) + map.getIm(1));
        v.assignMinus(abelMapAtZero);
        // v.assignDivide(2);
    }

    @Override
    public int getNumGenerators() {
        return numGenerators / 2;
    }
}
