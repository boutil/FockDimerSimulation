package dimerSim;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import lattices.VisualizationZ2;
import lattices.Z2Lattice;



public class LatticeTestAztec {
    
    public static void main(String[] args) {
        MarkovSimZ2 sim;

        // SchottkyDimersQuad schottkyDimers = buildSchottkyDimers();
        // Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, 301, 301);

        Z2Lattice lattice = new Z2Lattice(1001, 1001);


        // // // System.out.println(Arrays.deepToString(lattice.faceWeights));

        sim = new MarkovSimZ2(lattice, false);

        // sim = loadSim("experimentExport/AztecDiamond1001UniformConverged.ser");
        // sim = loadSim("AztecDiamond301G1symmetric.ser");

        // sim.setLattice(lattice);

        // sim.simulate((int)5e6);
        sim.simulate(1000000);


        // saveSim(sim, "AztecDiamond301G1symmetrictimes03.ser");
        saveSim(sim, "experimentExport/AztecDiamond1001UniformConverged.ser");

        // System.out.println(Arrays.deepToString(lattice.flipFaceWeights));

        VisualizationZ2 vis = new VisualizationZ2(sim);

        // vis.visualizeAmoeba(schottkyDimers);

        // vis.visualizeSim(5);

        // vis.visualizeWeights();

        // vis.visualizeDiscreteAbelMap();

        // vis.visualizeThetaCrossRatio(-3.1691590245331893, 5.071383906313774);

        // vis.visualizeDimerConfiguration();

        try {
            vis.saveDimerConfPic("experimentExport/AzteDiamond1001UniformConverged.png");
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

    public static MarkovSimZ2 loadSim(String fileName) {
        MarkovSimZ2 sim = null;
        try {
            ObjectInputStream in;
            in = new ObjectInputStream(new FileInputStream(fileName));
            Z2Lattice lattice = (Z2Lattice) in.readObject();
            byte[][] faceStates = (byte[][]) in.readObject();
            boolean[][] insideBoundary = (boolean[][]) in.readObject();
            sim = new MarkovSimZ2(lattice, faceStates, insideBoundary);
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

    public static MarkovSimZ2 simAztecDiamond(int size, int numSteps) {
        
        // We can do this to construct a simple Lattice and test simulation on it.
        int N = size, M = size;
        Z2Lattice lattice = new Z2Lattice(N, M);

        MarkovSimZ2 sim = new MarkovSimZ2(lattice, false);
        // sim.initializeAztecDiamond();

        long preSim = System.currentTimeMillis();
        // sim.simulate(numSteps);
        System.out.println("Simulation took " + (System.currentTimeMillis() - preSim)/1000 + " seconds");

        return sim;
    }

}
