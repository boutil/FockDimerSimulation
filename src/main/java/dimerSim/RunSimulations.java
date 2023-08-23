package dimerSim;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersDoubleCoverUnitary;
import de.jtem.riemann.schottky.SchottkyDimersUnitary;

public class RunSimulations {

    public static void main(String[] args) {
        // Create SchottkyDimers
        String baseFolder = "experimentExport/Hexagon/";

        // String experimentName = "01_starfish/";
        // int size = 300;
        // int numSteps = 500000;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "02_largeBubble/";
        // int size = 300;
        // int numSteps = 500000;
        // int interval = 0;
        // int numTimes = 10;
        
        String experimentName = "03_convergenceToUniform/";
        int size = 300;
        int numSteps = 500000;
        int interval = 5000;
        int numTimes = 1;
        
        
        String folderName = baseFolder + experimentName;
        
        // SchottkyDimers schottkyDimers = createSchottky();
        // SimulationManager man = new SimulationManager(schottkyDimers, folderName);

        continueSimulation(folderName, size, numSteps, interval, numTimes);
    }



    public static void continueSimulation(String folder, int size, int numSteps, int interval, int numTimes) {
        SimulationManager man = new SimulationManager(folder);

        man.setSavePictureInterval(interval);
        for (int i = 0; i < numTimes; i++) {
            man.simulateAndSave(numSteps, size);
        }
    }

    public static SchottkyDimers createSchottky(){
        Complex A = new Complex(Math.cos(0.1), Math.sin(0.1)).times(0.5);
        double[][] angles = {{0}, {Math.PI / 3}, {2 * Math.PI / 3}};

        Complex B = A.invert().conjugate();
        double[] schottkyParams = {A.re, A.im, B.re, B.im, 0.0000000001, 0};

        double[] sizes = {1, 1, 1};
        // double[] sizes = {1, 0.8, 0.6};

        double[] boundaryResidues = {sizes[0], -sizes[1], sizes[2], -sizes[0], sizes[1], -sizes[2]};

        SchottkyDimersUnitary schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(schottkyParams), angles, boundaryResidues);

        SchottkyDimersDoubleCoverUnitary doubleCover = schottkyDimers.getSimpleDoubleCover();

        return doubleCover;
    }



}
