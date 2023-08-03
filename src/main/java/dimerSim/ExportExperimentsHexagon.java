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
import de.jtem.riemann.schottky.SchottkyDimersHex;
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
        
        // Unitary case G1
        Complex A = new Complex(0.5 * Math.cos(0.5), 0.5 * Math.sin(0.5));
        Complex B = A.invert().conjugate();
        double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.001, 0}};
        double[][][] angles = {{{0, Math.PI * 2 / 3 - 0.1}, {Math.PI * 2 / 3, Math.PI * 2 / 3 + 0.5}, {Math.PI * 4 / 3, Math.PI * 4 / 3}}};

        // G2LargeHole
        // double[][] schottkyParamsCol = {{-1, 1, -1, -1, 0.03, 0, 1, 1, 1, -1, 0.001, 0}};
        // double[][][] angles = {{{-3}, {-2}, {-1}, {1}, {2}, {3}}};

        int defaultNumSteps = 100000;
        int[] numSteps = new int[schottkyParamsCol.length];
        Arrays.fill(numSteps, defaultNumSteps);
        // int[] numSteps = {100000, 100000, 100000};

        String baseFolder = "experimentExport/Hexagon/";
        // String simToStartFrom = "experimentExport/Hexagon/hexagon500UniformConverged.ser";
        String simToStartFrom = "experimentExport/Hexagon/2023-08-03-16-20-52/sim0[500x500].ser";

        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm-ss");
        LocalDateTime now = LocalDateTime.now();
        baseFolder += dtf.format(now);

        System.out.println(new File(baseFolder).mkdirs());

        for (int i = 0; i < schottkyParamsCol.length; i++) {

            // SchottkyDimersHex schottkyDimers = new SchottkyDimersHex(new SchottkyData(schottkyParamsCol[i]), angles[i]);

            SchottkyDimersUnitary schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(schottkyParamsCol[i]), angles[i]);

            SchottkyDimersUnitary doubleCover = schottkyDimers.getDoubleCover();
            
            // sim = new MarkovSimZ2(lattice, false);

            sim = loadSim(simToStartFrom);
            
            HexLatticeFock lattice = new HexLatticeFock(schottkyDimers, sim.lattice.N, sim.lattice.M);

            sim.setLattice(lattice);

            System.out.println("lattice built");

            sim.simulate(numSteps[i]);
    
            Visualization vis = new Visualization(sim);

            try {
                String info = i + "[" + sim.lattice.N + "x" + sim.lattice.M + "]";
                saveSim(sim, baseFolder + "/sim" + info + ".ser");
                saveSchottky(schottkyDimers, baseFolder + "/schottky" + info + ".ser");
                vis.saveAmoebaPic(schottkyDimers, baseFolder + "/amoebaPic" + info + ".png");
                vis.saveAztecPic(doubleCover, baseFolder + "/aztecPic" + info + ".png");
                vis.saveDimerConfPic(doubleCover, baseFolder + "/dimerConf" + info + ".png");
                // vis.saveWeightsPic(baseFolder + "/weights" + info + ".png");
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

