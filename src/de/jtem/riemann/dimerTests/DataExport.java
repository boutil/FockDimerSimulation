package de.jtem.riemann.dimerTests;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.Schottky;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.blas.ComplexMatrix;

import  java.io.BufferedWriter;
import java.lang.StringBuilder;
import java.util.Arrays;
import java.io.FileWriter;

public class DataExport {

    // TODO handle infinity or inverted a,b as an interval going through infinity
    public static Complex[] getPointArrayOnRealLine(double a, double b, int numPoints){
        Complex[] res = new Complex[numPoints];
        // take values from equidistant sampling on circle C(0, 1) and mapping [1, i, -1] via Moebius to [b, infty, a].
        if (b < a) {
            double rotationAngle = Math.PI / numPoints;
            double theta = 0.0;
            for (int i = 0; i < numPoints; i++) {
                double phi_theta = Math.signum(theta - Math.PI / 2) * Math.sqrt((1 + Math.sin(theta)) / (1 - Math.sin(theta)));
                res[i] = new Complex((a-b) / 2 * phi_theta + ((a + b) / 2), 0);
                theta += rotationAngle;
            }
            // for (int i = numPoints/2; i < numPoints; i++) {
            //     res[i] = new Complex(b - Math.pow(numPoints - i, 2), 0);
            // }
        } else {
            for (int i = 0; i < numPoints; i++) {
                res[i] = new Complex(a + ((b - a) / numPoints * i), 0);
            }
        }
        return res;
    }

    public static Complex[] getPointArrayOnCircle(Complex center, double radius, int numPoints) {
        Complex[] res = new Complex[numPoints];
        Complex rotation = Complex.exp(new Complex(0, 2 * Math.PI / (numPoints - 1)));
        Complex currentRotation = new Complex(radius, 0);
        for (int i = 0; i < numPoints; i++) {
            res[i] = new Complex(center.plus(currentRotation));
            currentRotation.assignTimes(rotation);
        }
        return res;
    }

    public static void main(String[] args) {
        Complex A = new Complex(-0.2, 0.4);
        Complex B = new Complex(0.5, 0.5);
        Complex muA = new Complex(0.01, 0);
        Complex muB = new Complex(0.1, 0);
        int schottkyGenus = 2;
        SchottkyData schottkyData = new SchottkyData( schottkyGenus );
        schottkyData.setA( 0, A );
        schottkyData.setB( 0, A.conjugate() );
        schottkyData.setMu( 0, muA );
        schottkyData.setA( 1, B );
        schottkyData.setB( 1, B.conjugate() );
        schottkyData.setMu( 1, muB );

        // Why is this the right A,B initiation and not the one above????
        // schottkyData.setA( 0, A );
        // schottkyData.setB( 0, A.neg() );
        // schottkyData.setMu( 0, muA );
        // schottkyData.setA( 1, A.conjugate().neg() );
        // schottkyData.setB( 1, A.conjugate() );
        // schottkyData.setMu( 1, muA );
        System.out.println(schottkyData.isClassical());
        Schottky schottky = new Schottky(schottkyData);
        System.out.println(schottky.isDifferentialSeriesEvaluable());
        System.out.println(schottky.isIntegralSeriesEvaluable());
        // ComplexMatrix perMat = schottky.getPeriodMatrix();
        // System.out.println(perMat);
        // System.out.println(Arrays.toString(schottky.numOfElementsWithWordLength));

        double[] angles = {-1, 0, 1};
        SchottkyDimers dimers = new SchottkyDimers(schottkyData, angles);

        int length = dimers.getCenters().length;
        for (int i = 0; i < length; i++) {
            System.out.println("" + dimers.getCenterOfCircle(i) + ", " + dimers.getRadius(i));
        }

        // create array of points P to compute the amoebaMap for
        int numPointsPerSegment = 100;
        int numSegments = 3 + schottkyGenus;
        Complex[][] points = new Complex[numSegments][numPointsPerSegment];
        points[0] = getPointArrayOnRealLine(angles[0], angles[1], numPointsPerSegment);
        points[1] = getPointArrayOnRealLine(angles[1], angles[2], numPointsPerSegment);
        points[2] = getPointArrayOnRealLine(angles[2], angles[0], numPointsPerSegment);
        points[3] = getPointArrayOnCircle(dimers.getCenterOfCircle(0), dimers.getRadius(0) * 1.0001, numPointsPerSegment);
        points[4] = getPointArrayOnCircle(dimers.getCenterOfCircle(1), dimers.getRadius(1) * 1.0001, numPointsPerSegment);

        // Calculate AmoebaMaps
        Complex[][] pointsAmoebaMapped = new Complex[numSegments][numPointsPerSegment];
        for (int i = 0; i < points.length; i++){
            for (int j = 0; j < points[i].length; j++){
                try {
                    pointsAmoebaMapped[i][j] = dimers.amoebaMapHexGrid(points[i][j]);
                } catch (Exception e) {
                    System.out.println("While calculating pointsAmoebaMapped: " + e.getMessage());
                }
            }
        }

        // Export values as csv
        try {
            BufferedWriter br = new BufferedWriter(new FileWriter("exportedPoints_2.json"));
            // your code
            StringBuilder sb = new StringBuilder();
            
            sb.append("{ \"points\": [");
            for (int i = 0; i < points.length; i++){
                sb.append("[");
                for (int j = 0; j < points[i].length; j++){
                    sb.append("[" + points[i][j].re + ", " + points[i][j].im + "]");
                    if(j != points[i].length - 1){
                        sb.append(", ");
                    }
                }
                sb.append("]");
                if(i != points.length - 1) {
                    sb.append(",");
                }
            }
            sb.append("], \n");
            sb.append(" \"pointsAmoebaMapped\": [");
            for (int i = 0; i < pointsAmoebaMapped.length; i++){
                sb.append("[");
                for (int j = 0; j < pointsAmoebaMapped[i].length; j++){
                    sb.append("[" + pointsAmoebaMapped[i][j].re + ", " + pointsAmoebaMapped[i][j].im + "]");
                    if(j != pointsAmoebaMapped[i].length - 1){
                        sb.append(", ");
                    }
                }
                sb.append("]");
                if(i != pointsAmoebaMapped.length - 1) {
                    sb.append(",");
                }
            }
            sb.append("]}");
            br.write(sb.toString());
            br.close();
        } catch (Exception e) {
            System.out.println("exception :" + e.getMessage());
        }
    }
}
