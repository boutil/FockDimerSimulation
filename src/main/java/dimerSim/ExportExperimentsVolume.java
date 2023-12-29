package dimerSim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;

import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import lattices.Visualization;
import lattices.Z2Lattice;
import lattices.Z2LatticeFock;

public class ExportExperimentsVolume {
    public static void main(String[] args) {
        // Create a folder and save all simulation results in that folder for easy readout.
        MarkovSimZ2 sim;
        
        // G1LargeHoleSide
        // Ex1
        // double[][] schottkyParamsCol = {{0.9, 1, 0.9, -1, 0.02, 0}};
        // double[][][] angles = {{{-2.4}, {-0.4}, {0.4}, {2.4}}};
        // String baseFolder = "experimentExport/Volume/501/03VolEx1G1Side/mu02/";
        
        // G1LargeHoleCentral
        // Ex2
        double[][] schottkyParamsCol = {{0.15, 1, 0.15, -1, 0.04, 0}};
        double[][][] angles = {{{-2.4}, {-0.4}, {0.4}, {2.4}}};
        String baseFolder = "experimentExport/Volume/501/03VolEx2G1Central/mu04/";
        
        // G2LargeHoles
        // Ex3
        // double[][] schottkyParamsCol = {{1.2, 1.3, 1.2, -1.3, 0.08, 0, -0.4, 0.6, -0.4, -0.6, 0.03, 0}};
        // // double[][] schottkyParamsCol = {{1.2, 1.3, 1.2, -1.3, 0.05, 0, -0.4, 0.6, -0.4, -0.6, 0.05, 0}};
        // double[][][] angles = {{{-2.4}, {-0.4}, {0.4}, {2.4}}};
        // String baseFolder = "experimentExport/Volume/03VolEx3G2/";
        
        // double[][] schottkyParamsCol = {{0.19, 1, 0.19, -1, 0.3, 0}};
        // double[][][] angles = {{{-2.4}, {-0.4}, {0.4}, {2.4}}};
        // String baseFolder = "experimentExport/Volume/03VolEx4G1Large/";
        
        // G1LargeHoleSide
        // Ex5
        // double[][] schottkyParamsCol = {{0.9, 1, 0.9, -1, 0.07, 0}};
        // double[][][] angles = {{{-2.4}, {-0.4}, {0.4}, {2.4}}};
        // String baseFolder = "experimentExport/Volume/03VolEx5G1Side07/";


        int defaultNumSteps = 200;
        int[] numSteps = new int[schottkyParamsCol.length];
        Arrays.fill(numSteps, defaultNumSteps);
        // int[] numSteps = {100000, 100000, 100000};

        // String baseFolder = "experimentExport/Volume/";
        // String simToStartFrom = "experimentExport/Volume/501/03VolEx1G1Side/mu02/2023-12-28-16-03-46/sim0[501x501].ser";
        String simToStartFrom = "experimentExport/Volume/501/03VolEx2G1Central/mu04/2023-12-29-09-23-15/sim0[501x501].ser";
        // String simToStartFrom = "experimentExport/Volume/501/2023-12-27-15-13-45/sim0[501x501].ser";
        // String simToStartFrom = "experimentExport/Volume/501/03VolUniformConverged[501x501].ser";

        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm-ss");
        LocalDateTime now = LocalDateTime.now();
        baseFolder += dtf.format(now);

        new File(baseFolder).mkdirs();

        for (int i = 0; i < schottkyParamsCol.length; i++) {
            
            SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(new SchottkyData(schottkyParamsCol[i]), angles[i]);
            
            sim = loadSim(simToStartFrom);
            
            Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, sim.lattice.N, sim.lattice.M);
            
            // Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, 301, 301);
            // Z2Lattice lattice = new Z2Lattice(501, 501);
            
            // sim = new MarkovSimZ2(lattice, true);
            
            // sim.numThreads = 1;
            
            sim.setLattice(lattice);
            sim.simulate(numSteps[i], 0.3);
    
            Visualization vis = new Visualization(sim, schottkyDimers);
            
            try {
                String info = i + "[" + sim.lattice.N + "x" + sim.lattice.M + "]";
                saveSim(sim, baseFolder + "/sim" + info + ".ser");
                saveSchottky(schottkyDimers, baseFolder + "/schottky" + info + ".ser");
                vis.saveDimerConfPic(baseFolder + "/aztecPic" + info + ".png", false);
                // vis.saveDimerConfPic(baseFolder + "/dimerConf" + info + ".png", true);
                vis.saveAmoebaPic(schottkyDimers, baseFolder + "/amoebaPic" + info + ".png");
                vis.saveWeightsPic(baseFolder + "/weights" + info + ".png");
                // vis.visualizeSim(3, new File(baseFolder, "heightFn.png"));
                // vis.visualizeWeights();
                // vis.visualizeThetaCrossRatio(lattice.abelIncrementsRight[0][0], lattice.abelIncrementsTop[0][0]);
            } catch (IOException e) {
                // TODO: handle exception
            }

        }
    }

    public static void saveSchottky(SchottkyDimersQuad schottkyDimers, String fileName) {
        try {
            ObjectOutputStream out;
            out = new ObjectOutputStream(new FileOutputStream(fileName));
            out.writeObject(schottkyDimers.getUniformizationData());
            out.writeObject(schottkyDimers.angles);
            out.flush();
            out.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static void saveSim(MarkovSimZ2 sim, String fileName) {
        try {
            ObjectOutputStream out;
            out = new ObjectOutputStream(new FileOutputStream(fileName));
            out.writeObject(sim.lattice);
            out.writeObject(sim.faceStates);
            out.writeObject(sim.insideBoundary);
            out.flush();
            out.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static SchottkyDimersQuad loadSchottky(String fileName) {
        SchottkyDimersQuad schottkyDimers = null;
        try {
            ObjectInputStream in;
            in = new ObjectInputStream(new FileInputStream(fileName));
            double[] uniformizationData = (double[]) in.readObject();
            double[][] angles = (double[][]) in.readObject();
            schottkyDimers = new SchottkyDimersQuad(new SchottkyData(uniformizationData), angles);
            in.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return schottkyDimers;
    }

    public static MarkovSimZ2 loadSim(String fileName) {
        MarkovSimZ2 sim = null;
        try {
            ObjectInputStream in;
            in = new ObjectInputStream(new FileInputStream(fileName));
            Z2Lattice lattice = (Z2Lattice) in.readObject();
            byte[][] faceStates = (byte[][]) in.readObject();
            boolean[][] insideBoundary = (boolean[][]) in.readObject();
            sim = new MarkovSimZ2(lattice, faceStates, insideBoundary);
            sim.computeHeightFunctionFromFaceStates();
            sim.computeVolume();
            in.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return sim;
    }


    // public static SchottkyDimersQuad buildSchottkyDimers(double[] schottkyParams, double[][] angles) {
    //     // First build a SchottkyData.
    //     double a = Math.sqrt(2 / (1.5 + Math.sqrt(2)));
    //     // double a = 0.01;
    //     // double[] schottkyParams = new double[]{0, 1, 0, -1, 0.2, 0};
        
    //     // Choose angles in a way such that crossratio is 1:
    //     double x = a * (1 + Math.sqrt(2));
    //     double firstAngle = -x - (a/2);
    //     // double[] angles = {firstAngle, firstAngle + x, firstAngle + x + a, firstAngle + x + a + x};
        
    //     // double[] schottkyParams = new double[]{firstAngle - 0.5, 0, firstAngle + x + a + x + 0.5, 0, 0.00005, 0};
    //     // double[] schottkyParams = new double[]{-1, 1, -1, -1, 0.05, 0, 1, 1, 1, -1, 0.05, 0};
    //     SchottkyData schottkyData = new SchottkyData(schottkyParams);

    //     // A nice not axis aligned example:
    //     // double[] angles = {-2.414, -0.414, 0.414, 2.414};
    //     for (int i = 0; i < angles.length; i++) {
    //         angles[i] *= 0.3;
    //         angles[i] -= 0;
    //     }
    //     System.out.println("angles crossRatio is " + crossRatio(angles));
    //     // Create the corresponding schottkyDimers.
    //     SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(schottkyData, angles);
    //     return schottkyDimers;
    // }

    // public static double crossRatio(double[] vals){
    //     return ((vals[1] - vals[0]) * (vals[3] - vals[2]) / (vals[2] - vals[1])) / (vals[0] - vals[3]);
    // }
}

