package dimerSim;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersQuad;



public class LatticeTest {
    
    public static void main(String[] args) {
        MarkovSimZ2 sim;

        SchottkyDimersQuad schottkyDimers = buildSchottkyDimers();
        Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, 100, 100);

        // System.out.println(Arrays.deepToString(lattice.faceWeights));

        sim = new MarkovSimZ2(lattice);
        // sim.initializeFlatSquare();
        // sim.simulate(10000);

        sim = simAztecDiamond(100, 100000);
        // sim = loadSim("mySim.ser");

        // saveSim(sim, "mySim.ser");

        VisualizationZ2 vis = new VisualizationZ2(sim);

        vis.visualizeSim();

        // vis.visualizeWeights();
        
    }

    public static void saveSim(MarkovSimZ2 sim, String fileName) {
        try {
            ObjectOutputStream out;
            out = new ObjectOutputStream(new FileOutputStream(fileName));
            out.writeObject(sim);
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
            sim = (MarkovSimZ2) in.readObject();
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

        MarkovSimZ2 sim = new MarkovSimZ2(lattice);
        sim.initializeAztecDiamond();

        long preSim = System.currentTimeMillis();
        sim.simulate(numSteps);
        System.out.println("Simulation took " + (System.currentTimeMillis() - preSim)/1000 + " seconds");

        return sim;
    }

    // public static MarkovSim2 simFockLattice()

    public static SchottkyDimersQuad buildSchottkyDimers() {
        // First build a SchottkyData.
        int schottkyGenus = 1;
        Complex mu = new Complex(0.3, 0);
        Complex A = new Complex(0.5, 1);
        // Complex B = new Complex();
        Complex B = A.conjugate();
        double[] schottkyParams = new double[]{A.re, A.im, B.re, B.im, mu.re, mu.im};
        SchottkyData schottkyData = new SchottkyData(schottkyParams);
        // Choose some angles
        double[] angles = {-5, -3, 1, 2};

        // Create the corresponding schottkyDimers.
        SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(schottkyData, angles);
        return schottkyDimers;
    }

}
