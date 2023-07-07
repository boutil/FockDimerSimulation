package dimerSim;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersQuad;

public class ExportExperiments {
    public static void main(String[] args) {
        // Create a folder and save all simulation results in that folder for easy readout.
        MarkovSimZ2 sim;

        SchottkyDimersQuad schottkyDimers = buildSchottkyDimers();
        Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, 301, 301);

        // sim = new MarkovSimZ2(lattice, false);

        sim = loadSim("AztecDiamond301UniformConverged.ser");
        // sim = loadSim("AztecDiamond301G1symmetric.ser");

        sim.setLattice(lattice);

        // sim.simulate(500000);

        
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






    public static SchottkyDimersQuad buildSchottkyDimers() {
        // First build a SchottkyData.
        double a = Math.sqrt(2 / (1.5 + Math.sqrt(2)));
        // double a = 0.01;
        double[] schottkyParams = new double[]{0, 1, 0, -1, 0.2, 0};
        
        // Choose angles in a way such that crossratio is 1:
        double x = a * (1 + Math.sqrt(2));
        double firstAngle = -x - (a/2);
        double[] angles = {firstAngle, firstAngle + x, firstAngle + x + a, firstAngle + x + a + x};
        
        // double[] schottkyParams = new double[]{firstAngle - 0.5, 0, firstAngle + x + a + x + 0.5, 0, 0.00005, 0};
        // double[] schottkyParams = new double[]{-1, 1, -1, -1, 0.05, 0, 1, 1, 1, -1, 0.05, 0};
        SchottkyData schottkyData = new SchottkyData(schottkyParams);

        // A nice not axis aligned example:
        // double[] angles = {-2.414, -0.414, 0.414, 2.414};
        for (int i = 0; i < angles.length; i++) {
            angles[i] *= 0.3;
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

