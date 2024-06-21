package inUtil;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersUnitary;

public class hexagonalCaseCheck {
    public static void main(String[] args) {
        Complex A = new Complex(Math.cos(2 * Math.PI / 3 + 0.1), Math.sin(2 * Math.PI / 3 + 0.1)).times(0.5);
        double[][][] angles = {{{0, 0}, {Math.PI / 3 - 0.6, Math.PI / 3 + 0.6}, {2 * Math.PI / 3, 2 * Math.PI / 3}}};
        // Complex A = new Complex(Math.cos(2 * Math.PI / 3 + 0.1), Math.sin(2 * Math.PI / 3 + 0.1)).times(0.5);
        // double[][][] angles = {{{0}, {Math.PI / 3}, {2 * Math.PI / 3}}};

        Complex B = A.invert().conjugate();
        double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.2, 0}};

        double[] boundaryResidues = {1, -1, 1, -1, 1, -1};
        // double[] boundaryResidues = {sizes[0], sizes[1], sizes[2], -sizes[0], -sizes[1], -sizes[2]};
        // double[] boundaryResidues = {-sizes[0], sizes[1], sizes[2], sizes[0], -sizes[1], -sizes[2]};
        
        SchottkyDimersUnitary schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(schottkyParamsCol[0]), angles[0], boundaryResidues);

        ComplexMatrix BMatrix = schottkyDimers.getPeriodMatrix();
        Complex pointOnInnerOval = schottkyDimers.getCenterOfCircle(0).plus(schottkyDimers.getRadius(0));
        Complex pointOnOuterOval = pointOnInnerOval.invert().conjugate();
        ComplexVector innerVec = new ComplexVector(1);
        ComplexVector outerVec = new ComplexVector(1);
        schottkyDimers.abelMap(innerVec, pointOnInnerOval);
        schottkyDimers.abelMap(outerVec, pointOnOuterOval);
        ComplexVector res = outerVec.minus(innerVec);
        Complex diff = res.get(0).minus(BMatrix.get(0, 0));
        System.out.println(diff);
    }
}
