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

        // This gives holes in sim.
        // Complex A = new Complex(Math.cos(3 * Math.PI / 6 - 0.5), Math.sin(3 * Math.PI / 6 - 0.5)).times(0.5);
        // Starfish:
        // Complex A = new Complex(Math.cos(Math.PI / 3 + 0.01), Math.sin(Math.PI / 3 + 0.01)).times(0.5);
        // double d = 0.07;
        // double[][][] angles = {{{d, Math.PI / 3 - d}, {Math.PI / 3 + d, 2 * Math.PI / 3 - d}, {2 * Math.PI / 3 + d, 3 * Math.PI / 3 - d}}};

        // One bubble:
        Complex A = new Complex(Math.cos(2 * Math.PI / 3 + 0.1), Math.sin(2 * Math.PI / 3 + 0.1)).times(0.5);
        double[][][] angles = {{{0, 0}, {Math.PI / 3 - 0.6, Math.PI / 3 + 0.6}, {2 * Math.PI / 3, 2 * Math.PI / 3}}};
        // Complex A = new Complex(Math.cos(2 * Math.PI / 3 + 0.1), Math.sin(2 * Math.PI / 3 + 0.1)).times(0.5);
        // double[][][] angles = {{{0}, {Math.PI / 3}, {2 * Math.PI / 3}}};

        Complex B = A.invert().conjugate();
        double[][] schottkyParamsCol = {{A.re, A.im, B.re, B.im, 0.06, 0}};
        
        
        int defaultNumSteps = 3000;
        int[] numSteps = new int[schottkyParamsCol.length];
        Arrays.fill(numSteps, defaultNumSteps);
        // int[] numSteps = {100000, 100000, 100000};
        
        String baseFolder = "experimentExport/Hexagon/";
        // String simToStartFrom = "experimentExport/Hexagon/hexagon300UniformConverged.ser";
        // String simToStartFrom = "experimentExport/Hexagon/hexagon400UniformConverged1-08-06.ser";
        // String simToStartFrom = "experimentExport/Hexagon/2023-08-23-09-56-41/sim0[300x300].ser";
        String simToStartFrom = "experimentExport/Hexagon/2023-08-23-11-23-40/sim0[300x300].ser";
        // String simToStartFrom = "experimentExport/Hexagon/2023-08-23-07-23-56/sim0[400x400].ser";
        
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm-ss");
        LocalDateTime now = LocalDateTime.now();
        baseFolder += dtf.format(now);

        System.out.println(new File(baseFolder).mkdirs());
        
        for (int i = 0; i < schottkyParamsCol.length; i++) {
            
            // SchottkyDimersHex schottkyDimers = new SchottkyDimersHex(new SchottkyData(schottkyParamsCol[i]), angles[i]);
            double[] sizes = {1, 1, 1};
            // double[] sizes = {1, 0.8, 0.6};

            double[] boundaryResidues = {sizes[0], -sizes[1], sizes[2], -sizes[0], sizes[1], -sizes[2]};
            // double[] boundaryResidues = {sizes[0], sizes[1], sizes[2], -sizes[0], -sizes[1], -sizes[2]};
            // double[] boundaryResidues = {-sizes[0], sizes[1], sizes[2], sizes[0], -sizes[1], -sizes[2]};
            
            SchottkyDimersUnitary schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(schottkyParamsCol[i]), angles[i], boundaryResidues);
            
            SchottkyDimersDoubleCoverUnitary doubleCover = schottkyDimers.getSimpleDoubleCover();
            
            // HexLatticeFock lattice = new HexLatticeFock(doubleCover, 400, 400);
            HexLattice lattice = new HexLattice(300, 300);
            sim = new MarkovSimHex(lattice, sizes);
            
            // sim = loadSim(simToStartFrom);
            
            // HexLatticeFock lattice = new HexLatticeFock(doubleCover, sim.lattice.N, sim.lattice.M);
            // HexLattice lattice = new HexLattice(sim.lattice.N, sim.lattice.M);

            // sim.setLattice(lattice);
            


            System.out.println("lattice built");

            // sim.simulate(numSteps[i]);
            sim.simulateGPU(numSteps[i]);
    
            Visualization vis = new Visualization(sim, doubleCover);

            try {
                String info = i + "[" + sim.lattice.N + "x" + sim.lattice.M + "]";
                saveSim(sim, baseFolder + "/sim" + info + ".ser");
                saveSchottky(doubleCover, baseFolder + "/schottky" + info + ".ser");
                vis.saveDimerConfPic(baseFolder + "/dimerConf" + info + ".png");
                vis.saveWeightsPic(baseFolder + "/weights" + info + ".png");
                vis.saveDimerConfPic(doubleCover, baseFolder + "/dimerConfPred" + info + ".png");
                vis.saveAztecPic(doubleCover, baseFolder + "/aztecPic" + info + ".png");
                vis.saveAmoebaPic(doubleCover, baseFolder + "/amoebaPic" + info + ".png");
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

