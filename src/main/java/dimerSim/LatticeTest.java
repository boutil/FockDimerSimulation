package dimerSim;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersQuad;



public class LatticeTest {
    
    public static void main(String[] args) {
        MarkovSimZ2 sim;

        SchottkyDimersQuad schottkyDimers = buildSchottkyDimers();
        Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, 200, 200);

        // Z2Lattice lattice = new Z2Lattice(400, 400);


        // // System.out.println(Arrays.deepToString(lattice.faceWeights));

        sim = new MarkovSimZ2(lattice);

        // sim = loadSim("mySimv2.ser");

        // sim.simulate(100000, 0.25);


        // saveSim(sim, "mySimCentralCase.ser");

        // System.out.println(Arrays.deepToString(lattice.flipFaceWeights));

        VisualizationZ2 vis = new VisualizationZ2(sim);

        // vis.visualizeSim(5);

        vis.visualizeWeights();

        // vis.visualizeDiscreteAbelMap();

        // vis.visualizeThetaCrossRatio(-3.1691590245331893, 5.071383906313774);
        
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
        // sim.initializeAztecDiamond();

        long preSim = System.currentTimeMillis();
        // sim.simulate(numSteps);
        System.out.println("Simulation took " + (System.currentTimeMillis() - preSim)/1000 + " seconds");

        return sim;
    }

    // public static MarkovSim2 simFockLattice()

    public static SchottkyDimersQuad buildSchottkyDimers() {
        // First build a SchottkyData.
        int schottkyGenus = 1;
        double a = Math.sqrt(2 / (1.5 + Math.sqrt(2)));
        // double a = 0.01;
        // double[] schottkyParams = new double[]{0, 1, 0, -1, 0.01, 0};
        
        // Choose angles in a way such that crossratio is 1:
        double x = a * (1 + Math.sqrt(2));
        double firstAngle = -x - (a/2);
        double[] angles = {firstAngle, firstAngle + x, firstAngle + x + a, firstAngle + x + a + x};
        
        // double[] schottkyParams = new double[]{firstAngle - 0.5, 0, firstAngle + x + a + x + 0.5, 0, 0.00005, 0};
        double[] schottkyParams = new double[]{-1, 1, -1, -1, 0.05, 0, 1, 1, 1, -1, 0.003, 0};
        SchottkyData schottkyData = new SchottkyData(schottkyParams);

        // A nice not axis aligned example:
        // double[] angles = {0, 0.1, 0.16, 0.61};
        for (int i = 0; i < angles.length; i++) {
            angles[i] *= 0.1;
            angles[i] -= 0;
        }
        System.out.println("angles crossRatio is " + crossRatio(angles));
        // Create the corresponding schottkyDimers.
        SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(schottkyData, angles);
        return schottkyDimers;
    }

    public static double crossRatio(double[] vals){
        return ((vals[1] - vals[0]) * (vals[3] - vals[2]) / (vals[2] - vals[1])) / (vals[0] - vals[3]);
    }

}
