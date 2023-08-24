package dimerSim;

import org.checkerframework.checker.units.qual.A;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersDoubleCoverUnitary;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import de.jtem.riemann.schottky.SchottkyDimersUnitary;

public class RunSimulations {

    public static void main(String[] args) {
        // Create SchottkyDimers
        String baseFolder = "experimentExport/";

        // String experimentName = "Hexagon/01_starfish/";
        // int size = 300;
        // int numSteps = 50;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Hexagon/02_largeBubble/";
        // int size = 300;
        // int numSteps = 500000;
        // int interval = 0;
        // int numTimes = 10;
        
        // String experimentName = "Hexagon/03_convergenceToUniform/";
        // int size = 300;
        // int numSteps = 60000;
        // int interval = 3000;
        // int numTimes = 1;
        
        // String experimentName = "Hexagon/04_starfishG0/";
        // int size = 300;
        // int numSteps = 500000;
        // int interval = 10000;
        // int numTimes = 1;
        
        
        // String experimentName = "Aztec/01_LargeHoleG1/";
        // int size = 501;
        // int numSteps = 5000000;
        // int interval = 250000;
        // int numTimes = 1;
        
        String experimentName = "Aztec/02_LargeHole2AnglesElongated/";
        int size = 100;
        int numSteps = 100;
        int interval = 0;
        int numTimes = 1;
        
        String folderName = baseFolder + experimentName;
        
        // SchottkyDimers schottkyDimers = createSchottkyHex();
        // SchottkyDimers schottkyDimers = createSchottkyQuad();


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

    public static SchottkyDimers createSchottkyHex(){
        // Starfish:
        Complex A = new Complex(Math.cos(Math.PI / 3 + 0.01), Math.sin(Math.PI / 3 + 0.01)).times(0.5);
        double d = 0.07;
        double[][] angles = {{d, Math.PI / 3 - d}, {Math.PI / 3 + d, 2 * Math.PI / 3 - d}, {2 * Math.PI / 3 + d, 3 * Math.PI / 3 - d}};

        Complex B = A.invert().conjugate();
        double[] schottkyParams = {A.re, A.im, B.re, B.im, 0.06, 0};

        double[] sizes = {1, 1, 1};
        // double[] sizes = {1, 0.8, 0.6};

        double[] boundaryResidues = {sizes[0], -sizes[1], sizes[2], -sizes[0], sizes[1], -sizes[2]};

        SchottkyDimersUnitary schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(schottkyParams), angles, boundaryResidues);

        SchottkyDimersDoubleCoverUnitary doubleCover = schottkyDimers.getSimpleDoubleCover();

        return doubleCover;
    }


    public static SchottkyDimers createSchottkyQuad(){
        // G1LargeHole
        double[] schottkyParamsCol = {0.9, 1, 0.9, -1, 0.15, 0};
        double[][] angles = {{-2.4}, {-0.4}, {0.4}, {2.4}};

        // G1LargeHole2Angles
        // double[] schottkyParamsCol = {0.9, 1, 0.9, -1, 0.15, 0};
        // double[][] angles = {{-2.4, -0.5}, {-0.4, -0.4}, {0.5, 1.3}, {1.4, 1.4}};

        SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(new SchottkyData(schottkyParamsCol), angles);

        return schottkyDimers;
    }


}
