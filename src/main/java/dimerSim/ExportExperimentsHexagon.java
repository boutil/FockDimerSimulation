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

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersDoubleCoverUnitary;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import de.jtem.riemann.schottky.SchottkyDimersUnitary;
import lattices.HexLattice;
import lattices.HexLatticeFock;
import lattices.Visualization;

public class ExportExperimentsHexagon {
    public static void main(String[] args) {
        // Create a folder and save all simulation results in that folder for easy readout.
        MarkovSimHex sim;

        // double[][] schottkyParamsCol = {{0.9, 1, 0.9, -1, 0.001, 0}};
        // double[][][] angles = {{{-3, -2}, {0.4, 0.4}, {2, 3}}};
        
        // // // Unitary case G1
        // Complex A = new Complex(Math.cos(Math.PI), Math.sin(Math.PI)).times(0.2);
        // Complex B = A.invert().conjugate();
        // double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.01, 0}};
        // // double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.05, 0}};
        // // double[][][] angles = {{{Math.PI * 4 / 3, Math.PI * 4 / 3}, {Math.PI * 2 - 0.5, Math.PI * 2 / 3 - 0.3}, {Math.PI * 2 / 3, Math.PI * 2 / 3 + 0.5}}};
        // double[][][] angles = {{{4 * Math.PI / 3}, {0}, {2 * Math.PI / 3}}};
        // // double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.00000000001, 0}};
        // // double[][][] angles = {{{Math.PI + 0.3, Math.PI + 0.3}, {Math.PI + 0.5, Math.PI - 0.5}, {Math.PI - 0.3, Math.PI - 0.3}}};
        
        // Complex A = new Complex(0.3 * Math.cos(0), 0.3 * Math.sin(0));
        // Complex B = A.invert().conjugate();
        // // double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.01, 0}};
        // double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.00000001, 0}};
        // double[][][] angles = {{{0, 0.5}, {Math.PI * 2 / 3, Math.PI * 2 / 3}, {Math.PI * 4 / 3, Math.PI * 4 / 3}}};

        // Complex A = new Complex(0.5 * Math.cos(0.5), 0.5 * Math.sin(0.5));
        // Complex B = A.invert().conjugate();
        // // double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.01, 0}};
        // double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.00000001, 0}};
        // double[][][] angles = {{{0, Math.PI * 2 / 3 - 0.08}, {Math.PI * 2 / 3, Math.PI * 2 / 3 + 0.5}, {Math.PI * 4 / 3, Math.PI * 4 / 3}}};

        // This gives holes in sim.
        Complex A = new Complex(Math.cos(Math.PI / 3 - 0.1), Math.sin(Math.PI / 3 - 0.1)).times(0.5);
        double[][][] angles = {{{0, 0}, {Math.PI / 3, Math.PI / 3}, {2 * Math.PI / 3 - 0.7, 2 * Math.PI / 3 + 0.7}}};
        // Complex A = new Complex(Math.cos(Math.PI / 6), Math.sin(Math.PI / 6)).times(0.6);
        // double[][][] angles = {{{0}, {Math.PI / 3}, {2 * Math.PI / 3}}};

        Complex B = A.invert().conjugate();
        double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.03, 0}};
        
        
        int defaultNumSteps = 500000;
        int[] numSteps = new int[schottkyParamsCol.length];
        Arrays.fill(numSteps, defaultNumSteps);
        // int[] numSteps = {100000, 100000, 100000};
        
        String baseFolder = "experimentExport/Hexagon/";
        String simToStartFrom = "experimentExport/Hexagon/hexagon300UniformConverged.ser";
        // String simToStartFrom = "experimentExport/Hexagon/2023-08-18-11-54-37/sim0[500x500].ser";
        // String simToStartFrom = "experimentExport/Hexagon/2023-08-18-13-34-22/sim0[300x300].ser";
        
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm-ss");
        LocalDateTime now = LocalDateTime.now();
        baseFolder += dtf.format(now);

        System.out.println(new File(baseFolder).mkdirs());

        for (int i = 0; i < schottkyParamsCol.length; i++) {

            // SchottkyDimersHex schottkyDimers = new SchottkyDimersHex(new SchottkyData(schottkyParamsCol[i]), angles[i]);
            // double[] sizes = {1, 1, 1};
            double[] sizes = {1, 0.8, 0.6};

            double[] boundaryResidues = {sizes[0], -sizes[1], sizes[2], -sizes[0], sizes[1], -sizes[2]};
            
            SchottkyDimersUnitary schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(schottkyParamsCol[i]), angles[i], boundaryResidues);

            SchottkyDimersDoubleCoverUnitary doubleCover = schottkyDimers.getSimpleDoubleCover();
            
            // HexLatticeFock lattice = new HexLatticeFock(doubleCover, 400, 400);
            // sim = new MarkovSimHex(lattice, sizes);
            
            sim = loadSim(simToStartFrom);
            
            HexLatticeFock lattice = new HexLatticeFock(doubleCover, sim.lattice.N, sim.lattice.M);
            // HexLattice lattice = new HexLattice(sim.lattice.N, sim.lattice.M);
            sim.setLattice(lattice);
            


            System.out.println("lattice built");

            sim.simulate(numSteps[i]);
    
            Visualization vis = new Visualization(sim);

            try {
                String info = i + "[" + sim.lattice.N + "x" + sim.lattice.M + "]";
                saveSim(sim, baseFolder + "/sim" + info + ".ser");
                saveSchottky(doubleCover, baseFolder + "/schottky" + info + ".ser");
                vis.saveDimerConfPic(baseFolder + "/dimerConf" + info + ".png");
                vis.saveWeightsPic(baseFolder + "/weights" + info + ".png");
                vis.saveDimerConfPic(doubleCover, baseFolder + "/dimerConfPred" + info + ".png");
                vis.saveAmoebaPic(doubleCover, baseFolder + "/amoebaPic" + info + ".png");
                // vis.saveAztecPic(doubleCover, baseFolder + "/aztecPic" + info + ".png");
            } catch (IOException e) {
                // TODO: handle exception
            }

        }
    }

    public static void saveSchottky(SchottkyDimers schottkyDimers, String fileName) {
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

    public static void saveSim(MarkovSimHex sim, String fileName) {
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

    public static MarkovSimHex loadSim(String fileName) {
        MarkovSimHex sim = null;
        try {
            ObjectInputStream in;
            in = new ObjectInputStream(new FileInputStream(fileName));
            HexLattice lattice = (HexLattice) in.readObject();
            byte[][] faceStates = (byte[][]) in.readObject();
            boolean[][] insideBoundary = (boolean[][]) in.readObject();
            sim = new MarkovSimHex(lattice, faceStates, insideBoundary);
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
}

