package dimerSim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersQuad;

public class ExportExperiments {
    public static void main(String[] args) {
        // Create a folder and save all simulation results in that folder for easy readout.
        MarkovSimZ2 sim;

        new File("asd").mkdirs();

        double[][] schottkyParamsCol = {new double[]{0, 1, 0, -1, 0.05, 0}, new double[]{0, 1, 0, -1, 0.2, 0}, new double[]{0.9, 1, 0.9, -1, 0.2, 0}};
        double[][] angles = {new double[]{-2.4, -0.4, 0.4, 2.4}, new double[]{-2.4, -0.4, 0.4, 2.4}, new double[]{-2.4, -0.4, 0.4, 2.4}};
        // int[] numSteps = {100000, 100000, 100000};
        int[] numSteps = {1000, 1000, 1000};

        String baseFolder = "experimentExport/";
        String simToStartFrom = "experimentExport/AztecDiamond301UniformConverged.ser";

        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm-ss");
        LocalDateTime now = LocalDateTime.now();
        baseFolder += dtf.format(now);

        System.out.println(new File(baseFolder).mkdirs());

        for (int i = 0; i < schottkyParamsCol.length; i++) {

            SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(new SchottkyData(schottkyParamsCol[i]), angles[i]);
            Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, 301, 301);
    
            // sim = new MarkovSimZ2(lattice, false);
    
            sim = loadSim(simToStartFrom);
            sim.setLattice(lattice);
            sim.simulate(numSteps[i]);
    
            VisualizationZ2 vis = new VisualizationZ2(sim);

            try {
                String info = i + "[" + lattice.N + "x" + lattice.M + "]";
                saveSim(sim, baseFolder + "/sim" + info + ".ser");
                vis.saveAmoebaPic(schottkyDimers, baseFolder + "/amoebaPic" + info + ".png");
                vis.saveDimerConfPic(baseFolder + "/dimerConf" + info + ".png");
                vis.saveWeightsPic(baseFolder + "/weights" + info + ".png");
            } catch (IOException e) {
                // TODO: handle exception
            }

        }
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






    public static SchottkyDimersQuad buildSchottkyDimers(double[] schottkyParams, double[] angles) {
        // First build a SchottkyData.
        double a = Math.sqrt(2 / (1.5 + Math.sqrt(2)));
        // double a = 0.01;
        // double[] schottkyParams = new double[]{0, 1, 0, -1, 0.2, 0};
        
        // Choose angles in a way such that crossratio is 1:
        double x = a * (1 + Math.sqrt(2));
        double firstAngle = -x - (a/2);
        // double[] angles = {firstAngle, firstAngle + x, firstAngle + x + a, firstAngle + x + a + x};
        
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

