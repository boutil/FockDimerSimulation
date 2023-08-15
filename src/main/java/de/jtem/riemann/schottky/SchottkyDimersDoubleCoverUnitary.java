package de.jtem.riemann.schottky;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;

public class SchottkyDimersDoubleCoverUnitary extends SchottkyDimersUnitary{

    // angles here interpreted as the arg of unitary complex angles.

    private ComplexVector abelMapAtZero;
    
    public SchottkyDimersDoubleCoverUnitary(SchottkyData data, double[][] angles, double[] boundaryResidues) {
        super(data, angles);

        Complex P0 = chooseP0();

        amoebaMap = new AmoebaMapHex(this, P0, boundaryResidues);

        ComplexVector abelMap0 = new ComplexVector(2, 0, 0);
        abelMap(abelMap0, new Complex(1,0));
        abelMapAtZero = new ComplexVector(1, abelMap0.getRe(0) + abelMap0.getRe(1),  abelMap0.getIm(0) + abelMap0.getIm(1));

        // It is assumed that the points in data and angles are doubled in z -> -z.
        adjustAngles();
        checkMCurveProp();
        checkSlopeAtWinding();
    }

    public SchottkyDimersDoubleCoverUnitary(SchottkyData data, double[][] angles) {
        this(data, angles, new double[] {1, -1, 1, -1, 1, -1});
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
        Complex slopeP0 = amoebaMap.getSlope(chooseP0(), acc);
        Complex boundaryDiff0 = amoebaMap.getDifferentials(new Complex(), acc)[2];
        Complex pointBetweenGammaAlpha = anglesC[2][anglesC[2].length - 1].plus(anglesC[3][0]);
        pointBetweenGammaAlpha.assignDivide(pointBetweenGammaAlpha.abs() / 0.99999);
        Complex pointBetweenBetaGamma = anglesC[1][anglesC[1].length - 1].plus(anglesC[2][0]);
        pointBetweenBetaGamma.assignDivide(pointBetweenBetaGamma.abs() / 0.99999);
        Complex slopeGammaAlpha = amoebaMap.getSlope(pointBetweenGammaAlpha, acc);
        Complex slopeBetaGamma = amoebaMap.getSlope(pointBetweenBetaGamma, acc);
        Complex correctSlope = slopeP0.plus(slopeGammaAlpha).plus(slopeBetaGamma).divide(3);
        System.out.println("uniformization data is: " + Arrays.toString(uniformizationData));
        System.out.println("angles are: " + Arrays.deepToString(angles));
        System.out.println("slope at P0 is " + slopeP0);
        System.out.println("slope at gammaAlpha is " + slopeGammaAlpha);
        System.out.println("slope at betaGamma is " + slopeBetaGamma);
        System.out.println("slope at inner oval is: " + amoebaMap.getSlope(getA(0).plus(getRadius(0)), acc));
        System.out.println("boundaryDifferential at 0 is: " + boundaryDiff0);
        System.out.println("slope at 0 should be: " + correctSlope);
        System.out.println("slope at 0 is " + amoebaMap.getSlope(new Complex(), acc));
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

    private void adjustAngles() {
        // adjusts angles so that dXiBoundary is 0 at 0.
        AmoebaMapHex amoebaMap = (AmoebaMapHex)this.amoebaMap;
        double eps = 0.01;
        double stepSize = 0.001;
        int maxSteps = 10000;
        int currentStep = 0;
        Complex[][] angleBdryDers = amoebaMap.calculateBoundaryDerivativesByAngles();
        while(amoebaMap.dXiBoundary.abs() > eps) {
            for (int i = 0; i < angleBdryDers.length; i++) {
                for (int j = 0; j < angleBdryDers[i].length; j++) {
                    // project the derivative onto -dXiBoundary to see by how much to move.
                    double projectionLength = amoebaMap.dXiBoundary.times(angleBdryDers[i][j].conjugate()).re / amoebaMap.dXiBoundary.abs();
                    angles[i][j] -= projectionLength * stepSize;
                    angles[i][j] = doubleMod(angles[i][j], 2 * Math.PI);
                    angles[i + 3][j] = doubleMod(angles[i][j] + Math.PI, 2 * Math.PI);
                }
            }
            angleBdryDers = amoebaMap.calculateBoundaryDerivativesByAngles();
            System.out.println(amoebaMap.dXiBoundary.abs());
        }
    }

    @Override
    public int getNumGenerators() {
        return numGenerators / 2;
    }
}
