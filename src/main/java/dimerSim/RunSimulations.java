package dimerSim;

import java.io.File;

import org.checkerframework.checker.units.qual.A;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersDoubleCoverUnitary;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import de.jtem.riemann.schottky.SchottkyDimersQuadUnitary;
import de.jtem.riemann.schottky.SchottkyDimersUnitary;

public class RunSimulations {

    public static void main(String[] args) {
        // Create SchottkyDimers
        String baseFolder = "experimentExport/";

        // String experimentName = "Hexagon/01_starfish/";
        // int size = 600;
        // int numSteps = (int) 1e7;
        // int interval = 0;
        // int numTimes = 10;
        
        // String experimentName = "Hexagon/02_growingLargeBubble/mu005/";
        // int size = 600;
        // int numSteps = (int) 1e7;
        // int interval = 0;
        // int numTimes = 10;
        
        // String experimentName = "Hexagon/03_convergenceToUniform/";
        // int size = 600;
        // int numSteps = (int)1e6;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Hexagon/04_growingStarfishHoles/mu06/";
        // int size = 300;
        // int numSteps = (int)5e6;
        // int interval = 0;
        // int numTimes = 10;

        // --------------------------------------------------

        // String experimentName = "Aztec/01_LargeHoleG1/";
        // int size = 604;
        // int numSteps = (int) 1e7;
        // int interval = (int) 0;
        // int numTimes = 10;
        
        // String experimentName = "Aztec/02_LargeHole2AnglesElongated/";
        // int size = 1000;
        // int numSteps = (int)1e6;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/03_uniform/";
        // String experimentName = "Aztec/05_growing1LargeHole/mu004/";
        // int size = 500;
        // int numSteps = (int) 1e7;
        // int interval = 0;
        // int numTimes = 10;
        
        // String experimentName = "Aztec/06_growing2LargeHolesFromUnitary/mu01/";
        // int size = 600;
        // int numSteps = (int) 1e7;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/07_growing2LargeHolesFromUnitarySymmetric/mu01/";
        // int size = 600;
        // int numSteps = (int) 1e7;
        // int interval = 0;
        // int numTimes = 1;
        
        String experimentName = "Aztec/08_growing2LargeHolesFromUnitarySymmetric1Angle/mu02/";
        int size = 600;
        int numSteps = (int) 1e7;
        int interval = 0;
        int numTimes = 1;

        System.out.println("Running experiment " + experimentName + " of size " + size);
        
        String folderName = baseFolder + experimentName;
        
        // SchottkyDimers schottkyDimers = createSchottkyQuad();
        // SchottkyDimers schottkyDimers = createSchottkyHex();


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


    public static void growHoles(String baseFolder, int size, int numSteps, int interval, int numTimes, double[] mus) {
        for (int i = 1; i < mus.length; i++) {
            String muFolder = new File(baseFolder, String.format("%0.5f", mus[i])).toString();

        }
    }


    public static SchottkyDimers createSchottkyHex(){
        // Starfish:
        // Complex A = new Complex(Math.cos(Math.PI / 3 + 0.01), Math.sin(Math.PI / 3 + 0.01)).times(0.5);
        // double d = 0.07;
        // double[][] angles = {{d, Math.PI / 3 - d}, {Math.PI / 3 + d, 2 * Math.PI / 3 - d}, {2 * Math.PI / 3 + d, 3 * Math.PI / 3 - d}};

        // One bubble with angles
        Complex A = new Complex(Math.cos(2 * Math.PI / 3 + 0.1), Math.sin(2 * Math.PI / 3 + 0.1)).times(0.5);
        double[][] angles = {{0, 0}, {Math.PI / 3 - 0.6, Math.PI / 3 + 0.6}, {2 * Math.PI / 3, 2 * Math.PI / 3}};

        Complex B = A.invert().conjugate();
        double[] schottkyParams = {A.re, A.im, B.re, B.im, 0.05, 0};

        double[] sizes = {1, 1, 1};
        // double[] sizes = {1, 0.8, 0.6};

        double[] boundaryResidues = {sizes[0], -sizes[1], sizes[2], -sizes[0], sizes[1], -sizes[2]};

        SchottkyDimersUnitary schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(schottkyParams), angles, boundaryResidues);

        SchottkyDimersDoubleCoverUnitary doubleCover = schottkyDimers.getSimpleDoubleCover();

        return doubleCover;
    }


    public static SchottkyDimers createSchottkyQuad(){
        // G1LargeHole
        // double[] schottkyParamsCol = {0.9, 1, 0.9, -1, 0.08, 0};
        // double a = Math.sqrt(2 / (1.5 + Math.sqrt(2)));
        // double x = a * (1 + Math.sqrt(2));
        // double firstAngle = -x - (a/2);
        // double[][] angles = {{firstAngle}, {firstAngle + x}, {firstAngle + x + a}, {firstAngle + x + a + x}};

        // double[][] angles = {{-2.4}, {-0.4}, {0.6}, {2.4}};
        
        // G1LargeHole2Angles
        // double[] schottkyParamsCol = {0.9, 1, 0.9, -1, 0.15, 0};
        // double[][] angles = {{-2.4, -0.5}, {-0.4, -0.4}, {0.5, 1.3}, {1.4, 1.4}};
        
        // SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(new SchottkyData(schottkyParamsCol), angles);
        
        // // G1Unitary
        // double theta = Math.PI / 4;
        // Complex A = new Complex(Math.cos(theta), Math.sin(theta)).times(0.3);
        // Complex ARefl = A.invert().conjugate();
        // double[] schottkyParamsCol = {A.re, A.im, ARefl.re, ARefl.im, 0.00000004, 0};
        // double[][] angles = {{0}, {Math.PI / 2}, {Math.PI}, {3 * Math.PI / 2}};
        
        // G2Unitary
        double theta1 = Math.PI / 2 + 0.7;
        Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.3);
        Complex A1Refl = A1.invert().conjugate();
        double theta2 = -Math.PI / 4;
        Complex A2 = new Complex(Math.cos(theta2), Math.sin(theta2)).times(0.3);
        Complex A2Refl = A2.invert().conjugate();
        double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.02, 0, A2.re, A2.im, A2Refl.re, A2Refl.im, 0.02, 0};
        double[][] angles = {{0, 0}, {Math.PI / 2 - 0.9, Math.PI / 2 + 0.4}, {Math.PI, Math.PI}, {3 * Math.PI / 2, 3 * Math.PI / 2}};
        
        SchottkyDimersQuadUnitary schottkyDimers = new SchottkyDimersQuadUnitary(new SchottkyData(schottkyParamsCol), angles);
        return schottkyDimers;
    }


}
