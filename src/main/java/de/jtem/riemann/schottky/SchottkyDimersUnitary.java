package de.jtem.riemann.schottky;

import java.util.LinkedList;
import java.util.List;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

import de.jtem.mfc.field.Complex;

public class SchottkyDimersUnitary extends SchottkyDimers{

    // angles here interpreted as the arg of unitary complex angles.
    
    public SchottkyDimersUnitary(SchottkyData data, double[][] angles) {
        super(data, angles);

        // TODO: P0 should be chosen such that it's in the fundamental domain.
        Complex P0 = new Complex(0, 0);

        amoebaMap = new AmoebaMapHex(this, P0);

        checkMCurveProp();
    }

    private void checkMCurveProp() {
        double eps = 0.01;
        for (int i = 0; i < numGenerators; i++) {
            if (!getA(i).equals(getB(i).invert().conjugate(), eps)){
                System.out.println("not a unitary M-curve");
            }
        }
    }

    private double[] getAnglesInOrder() {
        double[] orderedAngles = new double[numAngles];
        for (int i = 0; i < angles.length; i++) {
            for (int j = 0; j < angles[i].length; j++) {
                orderedAngles[i + j] = angles[i][j];
            }
        }
        return orderedAngles;
    }

    @Override
    public Complex[] getPointArrayOnRealLine(double a, double b, int numPoints, double[] excludingInterval){
        List<Complex> res = new LinkedList<Complex>();
        double radius = 1.01;
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

    public SchottkyDimersUnitary getDoubleCover() {
        double[][] newAngles = new double[angles.length * 2][];
        for (int i = 0; i < angles.length; i++) {
            newAngles[i] = new double[angles[i].length];
            newAngles[i + angles.length] = new double[angles[i].length];
            for (int j = 0; j < angles[i].length; j++) {
                newAngles[i][j] = angles[i][j] / 2;
                newAngles[i + angles.length][j] = angles[i][j] / 2 + Math.PI;
            }
        }
        SchottkyData data = new SchottkyData(numGenerators * 2);
        for (int i = 0; i < numGenerators; i++) {
            double theta = getA(i).arg()/2;
            Complex rotation = new Complex(Math.cos(theta), Math.sin(theta));
            data.setA(i, getA(i).divide(rotation));
            data.setB(i, getB(i).divide(rotation));
            data.setMu(i, getMu(i));
            data.setA(i + numGenerators, getA(i).divide(rotation).neg());
            data.setB(i + numGenerators, getB(i).divide(rotation).neg());
            data.setMu(i + numGenerators, getMu(i));
        }
        return new SchottkyDimersUnitary(data, newAngles);
    }
}
