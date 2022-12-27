package de.jtem.riemann.dimerTests;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersHex;
import de.jtem.riemann.schottky.SchottkyDimersQuad;

public class AmoebaExportTest {
    public static void main(String[] args) {

        double a = Math.sqrt(2 / (1.5 + Math.sqrt(2)));
        // double[] schottkyParams = new double[]{-a, 1, -a, -1, 0.3, 0};
        double[] schottkyParams = new double[]{-a, 1, -a, -1, 0.15, 0, a, 1, a, -1, 0.01, 0};

        SchottkyData schottkyData = new SchottkyData(schottkyParams);
        
        // Choose some angles
        // double[] angles = {-0.7, -0.1, 0.4, 1};
        // Most symmetric angle setup.
        // double a = 0.05;
        double x = a * (1 + Math.sqrt(2));
        double firstAngle = -x - (a/2);
        double[] angles = {firstAngle, firstAngle + x, firstAngle + x + a, firstAngle + x + a + x};
        // Create the corresponding schottkyDimers.
        SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(schottkyData, angles);
        int numPointsPerSegment = 500;      
        Complex[][] points = schottkyDimers.parametrizeRealOvals(numPointsPerSegment);
        double[][][] amoebaPoints = extractAmoebaPointsFromSchottky(schottkyDimers, points);
        Map<String, Object> info = encodeAmoeba(schottkyDimers, amoebaPoints, points);

        
        List<Map<String, Object>> allAmoebas = new LinkedList<Map<String, Object>>();
        // for (int i = 0; i < params.length; i++) {
        //     // SchottkyData schottkyData = new SchottkyData( schottkyGenus );
        //     double[] schottkyParams = new double[]{params[i].A.re, params[i].A.im, params[i].B.re, params[i].B.im, params[i].mu.re, params[i].mu.im};
        //     SchottkyData schottkyData = new SchottkyData( schottkyParams );
        //     // schottkyData.setA( 0, params[i].A );
        //     // schottkyData.setB( 0, params[i].B );
        //     // schottkyData.setMu( 0, params[i].mu );
        //     SchottkyDimersHex dimers = new SchottkyDimersHex(schottkyData, params[i].angles);
        //     Complex[][] points = dimers.parametrizeRealOvals(numPointsPerSegment);
        //     double[][][] amoebaPoints = extractAmoebaPointsFromSchottky(dimers, points);
        //     Map<String, Object> info = encodeAmoeba(dimers, amoebaPoints, points);
        //     allAmoebas.add(info);
        // }

        allAmoebas.add(info);


        Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
        String allAmoebasJSON = gson.toJson(allAmoebas);

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter("exportedPoints_QuadAmoeba.json"));
            bw.write(allAmoebasJSON);
            bw.close();
        } catch (Exception e) {
            System.out.println("exception :" + e.getMessage());
        }
    }

    public static class SchottkyDimersData {
        public Complex A;
        public Complex B;
        public Complex mu;
        public double[] angles;
        public SchottkyDimersData(Complex A, Complex B, Complex mu, double[] angles) {
            this.A = A;
            this.B = B;
            this.mu = mu;
            this.angles = angles;
        }
    }

    public static SchottkyDimersData[] generateSchottkyParams(double[] mus, double[][] ABs, double[][] angles) {
        List<SchottkyDimersData> params = new LinkedList<SchottkyDimersData>();
        for(double mu : mus) {
            for(double[] AB : ABs) {
                for(double[] angle : angles){
                    // Does U1 parametrization
                    params.add(new SchottkyDimersData(new Complex(AB[0], 0), new Complex(AB[1], 0), new Complex(mu, 0), angle));
                }
            }
        }
        return params.toArray(SchottkyDimersData[]::new);
    }

    public static double[][][] extractAmoebaPointsFromSchottky(SchottkyDimers dimers, Complex[][] points) {
        // Not a very general form for now.
        int numSegments = dimers.angles.length + dimers.getNumGenerators();
        // Complex[][] points = dimers.parametrizeRealOvals(numPointsPerSegment);
        double[] angles = dimers.angles;
        double[][][] pointsAmoebaMapped = new double[numSegments][][];
        for (int i = 0; i < points.length; i++){
            List<Double[]> mappedP = new LinkedList<Double[]>();
            // double[][] mappedP = new double[points[i].length][2];
            for (int j = 0; j < points[i].length; j++){
                try {
                    mappedP.add(new Double[]{dimers.amoebaMap(points[i][j]).re,dimers.amoebaMap(points[i][j]).im});
                    // mappedP[j][0] = dimers.amoebaMapHexGrid(points[i][j]).re;
                    // mappedP[j][1] = dimers.amoebaMapHexGrid(points[i][j]).im;
                } catch (Exception e) {
                    // System.out.println("While calculating pointsAmoebaMapped: " + e.getMessage() + "P: " + points[i][j]);
                    // if (j != 0) { 
                        //     mappedP[j] = mappedP[j-1];
                    // }
                }
            }
            double[][] mappedArray = new double[mappedP.size()][2];
            int j = 0;
            for (Double[] d : mappedP) {
                mappedArray[j][0] = d[0];
                mappedArray[j][1] = d[1];
                j++;
            }
            // pointsAmoebaMapped[i] = mappedP;
            pointsAmoebaMapped[i] = mappedArray;
        }
        // values.add(points);
        // values.add(pointsAmoebaMapped);
        return pointsAmoebaMapped;
    }

    public static Map<String, Object> encodeAmoeba(SchottkyDimers dimers, double[][][] pointsAmoebaMapped, Complex[][] points) {
        Map<String, Object> info = new HashMap<String, Object>();
        info.put("pointsAmoebaMapped", pointsAmoebaMapped);
        info.put("angles", dimers.angles);
        info.put("uniformizationData", dimers.getUniformizationData());
        info.put("pointsOnOvals", points);
        return info;
    }
}
