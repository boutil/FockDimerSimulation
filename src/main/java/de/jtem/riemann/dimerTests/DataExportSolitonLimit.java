// package de.jtem.riemann.dimerTests;

// import java.io.BufferedWriter;
// import java.io.FileWriter;
// import java.util.HashMap;
// import java.util.LinkedList;
// import java.util.List;
// import java.util.Map;

// import com.google.gson.Gson;
// import com.google.gson.GsonBuilder;

// import de.jtem.mfc.field.Complex;
// import de.jtem.riemann.schottky.SchottkyData;
// import de.jtem.riemann.schottky.SchottkyDimersHex;

// public class DataExportSolitonLimit {
//     public static void main(String[] args) {
//         int schottkyGenus = 2;
//         double[][] mus = {{0.000000000000001, 0.01}, {0.01}};
//         // double[][] ABs = {{-0.51, -0.49}, {-0.75, -0.25}, {-0.99, -0.01}, {-3, -2}, {-3, -1.1}};
//         // double[][] ABs = {{-0.009999, -0.00001}, {-0.009999999, -0.000000001} };
//         double[][][] ABs = {{{-0.9, -0.1}, {0.3, 0.7}}, {{0.3, 0.7}}};

//         // double[][] ABs = {{0.2, 0.8}, {0.1, 0.9}, {0.01, 0.99}, {0.001, 0.999}, {0.0001, 0.9999}, {-0.8, -0.2}, {-0.9, -0.1}, {-0.99, -0.01}, {-0.999, -0.001}, {-0.9999, -0.0001}, {-0.99999, -0.00001}, {-0.999999, -0.000001}, {-0.9999999, -0.0000001}};
//         // double[][] angles = {{-1, 0, 1}, {-1, -0.5, 0}, {-1, -0.90, -0.8}, {-1, -0.99, -0.98}, {-1.9, -1.8, -1.7}};
//         // double[][] angles = {{-5, 5, 100}};
//         double[][] angles = {{-1, 0, 1}, {-1, 0, 1}};
//         SchottkyDimersData[] schottkyDimerDatas = generateSchottkyDimerDatas(ABs, mus, angles);
        
//         int numPointsPerSegment = 500;      
//         List<Map<String, Object>> allAmoebas = new LinkedList<Map<String, Object>>();
//         for (int i = 0; i < schottkyDimerDatas.length; i++) {
//             // SchottkyData schottkyData = new SchottkyData( schottkyGenus );
//             // schottkyData.setA( 0, params[i].A );
//             // schottkyData.setB( 0, params[i].B );
//             // schottkyData.setMu( 0, params[i].mu );
//             SchottkyDimersHex dimers = new SchottkyDimersHex(schottkyDimerDatas[i].schottkyData, schottkyDimerDatas[i].angles);
//             Complex[][] points = dimers.parametrizeRealOvals(numPointsPerSegment);
//             double[][][] amoebaPoints = extractAmoebaPointsFromSchottky(dimers, points);
//             Map<String, Object> info = encodeAmoeba(dimers, amoebaPoints, points);
//             allAmoebas.add(info);
//         }


//         Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
//         String allAmoebasJSON = gson.toJson(allAmoebas);

//         try {
//             BufferedWriter bw = new BufferedWriter(new FileWriter("exportedPoints_soliton_limit.json"));
//             bw.write(allAmoebasJSON);
//             bw.close();
//         } catch (Exception e) {
//             System.out.println("exception :" + e.getMessage());
//         }
//     }

//     public static class SchottkyDimersData {
//         public SchottkyData schottkyData;
//         public double[] angles;
//         // Assuming U1 Parametrization
//         public SchottkyDimersData(double[][] ABs, double[] mus, double[] angles) {
//             double[] schottkyParams = new double[mus.length * 6];
//             for (int i = 0, j = 0; i < mus.length; i++){
//                 schottkyParams[j++] = ABs[i][0];
//                 schottkyParams[j++] = 0;
//                 schottkyParams[j++] = ABs[i][1];
//                 schottkyParams[j++] = 0;
//                 schottkyParams[j++] = mus[i];
//                 schottkyParams[j++] = 0;
//             }
//             schottkyData = new SchottkyData(schottkyParams);
//             this.angles = angles;
//         }
//     }

//     public static SchottkyDimersData[] generateSchottkyDimerDatas(double[][][] ABs, double[][] mus, double[][] angles) {
//         // Here we assume that all arrays have the same length
//         List<SchottkyDimersData> params = new LinkedList<SchottkyDimersData>();
//         for (int i = 0; i < angles.length; i++){
//             params.add(new SchottkyDimersData(ABs[i], mus[i], angles[i]));
//         }
//         // for(double[] mu : mus) {
//         //     for(double[][] AB : ABs) {
//         //         for(double[] angle : angles){
//         //             // Does U1 parametrization
//         //             params.add(new SchottkyDimersData(AB, mu, angle));
//         //         }
//         //     }
//         // }
//         return params.toArray(SchottkyDimersData[]::new);
//     }

//     public static double[][][] extractAmoebaPointsFromSchottky(SchottkyDimersHex dimers, Complex[][] points) {
//         // Not a very general form for now.
//         int numSegments = 3 + dimers.getNumGenerators();
//         // Complex[][] points = dimers.parametrizeRealOvals(numPointsPerSegment);
//         double[] angles = dimers.angles;
//         double[][][] pointsAmoebaMapped = new double[numSegments][][];
//         for (int i = 0; i < points.length; i++){
//             List<Double[]> mappedP = new LinkedList<Double[]>();
//             // double[][] mappedP = new double[points[i].length][2];
//             for (int j = 0; j < points[i].length; j++){
//                 try {
//                     mappedP.add(new Double[]{dimers.amoebaMap(points[i][j]).re,dimers.amoebaMap(points[i][j]).im});
//                     // mappedP[j][0] = dimers.amoebaMapHexGrid(points[i][j]).re;
//                     // mappedP[j][1] = dimers.amoebaMapHexGrid(points[i][j]).im;
//                 } catch (Exception e) {
//                     // System.out.println("While calculating pointsAmoebaMapped: " + e.getMessage() + "P: " + points[i][j]);
//                     // if (j != 0) { 
//                         //     mappedP[j] = mappedP[j-1];
//                     // }
//                 }
//             }
//             double[][] mappedArray = new double[mappedP.size()][2];
//             int j = 0;
//             for (Double[] d : mappedP) {
//                 mappedArray[j][0] = d[0];
//                 mappedArray[j][1] = d[1];
//                 j++;
//             }
//             // pointsAmoebaMapped[i] = mappedP;
//             pointsAmoebaMapped[i] = mappedArray;
//         }
//         // values.add(points);
//         // values.add(pointsAmoebaMapped);
//         return pointsAmoebaMapped;
//     }

//     public static Map<String, Object> encodeAmoeba(SchottkyDimersHex dimers, double[][][] pointsAmoebaMapped, Complex[][] points) {
//         Map<String, Object> info = new HashMap<String, Object>();
//         info.put("pointsAmoebaMapped", pointsAmoebaMapped);
//         info.put("angles", dimers.angles);
//         info.put("uniformizationData", dimers.getUniformizationData());
//         info.put("pointsOnOvals", points);
//         return info;
//     }
// }


