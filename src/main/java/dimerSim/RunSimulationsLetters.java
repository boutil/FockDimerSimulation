package dimerSim;

import java.io.File;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersDoubleCoverUnitary;
import de.jtem.riemann.schottky.SchottkyDimersQuadUnitary;
import de.jtem.riemann.schottky.SchottkyDimersUnitary;

public class RunSimulationsLetters {

    public static void main(String[] args) {
        // Create SchottkyDimers
        String baseFolder = "experimentExport/";

        // String experimentName = "Hexagon/05_g0anglesOn2Sides/mu00000001/";
        // int size = 600;
        // int numSteps = (int)1e2;
        // int interval = 0;
        // int numTimes = 1;

        // --------------------------------------------------
        int size = 600;

        String experimentName = "AztecLetters/Letter_D";
        // String experimentName = "AztecLetters/Letter_I";
        // String experimentName = "AztecLetters/Letter_M";
        // String experimentName = "AztecLetters/Letter_E";
        // String experimentName = "AztecLetters/Letter_R";
        // String experimentName = "AztecLetters/Letter_S";


        System.out.println("Running experiment " + experimentName + " of size " + size);
        
        String folderName = baseFolder + experimentName;
        
        SchottkyDimers schottkyDimers = createSchottkyQuad();
        // // SchottkyDimers schottkyDimers = createSchottkyHex();

        SimulationManager man = new SimulationManager(schottkyDimers, folderName);

        // continueSimulation(folderName, size, numSteps, interval, numTimes);
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


    public static SchottkyDimers createSchottkyQuad(){
        // Letter_D
        double theta1 = - Math.PI / 4 - 0.6;
        Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.1);
        Complex A1Refl = A1.invert().conjugate();
        double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.1, 0};
        double[][] angles = {{0.5, 0.5}, {Math.PI / 2 + 0.3, Math.PI / 2 + 0.6}, {Math.PI, Math.PI}, {3 * Math.PI / 2, 3 * Math.PI / 2}};
        
        // // Letter_I
        // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.1);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000000001, 0};
        // double[][] angles = {{0}, {Math.PI / 2 - 1.2}, {Math.PI}, {3 * Math.PI / 2}};

        // // Letter_M
        // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.1);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000000001, 0};
        // double[][] angles = {{0, 0, 0, 0}, {Math.PI / 2 - 1.3, Math.PI / 2 - 1.3, Math.PI / 2 + 1.3, Math.PI / 2 + 1.3}, {Math.PI, Math.PI, Math.PI, Math.PI}, {3 * Math.PI / 2- 1.49, 3 * Math.PI / 2, 3 * Math.PI / 2, 3 * Math.PI / 2 + 1.49}};

        // // Letter_E
        // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.1);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000000001, 0};
        // double[][] angles = {{0.7, 0.7, 0.7}, {Math.PI / 2, Math.PI / 2, Math.PI / 2}, {Math.PI - 0.7, Math.PI - 0.7, Math.PI - 0.7}, {3 * Math.PI / 2 - 2, 3 * Math.PI / 2, 3 * Math.PI / 2 + 2}};
        
        // // Letter_R
        // double theta1 = Math.PI / 2 - 0.1;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.4);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.03, 0};
        // double[][] angles = {{0, 0, 0}, {Math.PI / 2 - 0.3, Math.PI / 2 - 0.3, Math.PI / 2 - 0.3}, {Math.PI - 1.3, Math.PI + 1.45, Math.PI + 1.45}, {3 * Math.PI / 2, 3 * Math.PI / 2 + 1.2, 3 * Math.PI / 2 + 1.2}};

        // // Letter_S
        // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.1);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000000001, 0};
        // double[][] angles = {{-0.1, -0.1, -0.1, 1.6, 1.6}, {Math.PI / 2 + 0.25, Math.PI / 2 + 0.25, Math.PI / 2 + 0.25, Math.PI / 2 + 0.25, Math.PI / 2 + 0.25}, {Math.PI - 0.5, Math.PI - 0.5, Math.PI - 0.5, Math.PI + 1.49, Math.PI + 1.49}, {3 * Math.PI / 2 + 0.25, 3 * Math.PI / 2 + 0.25, 3 * Math.PI / 2 + 0.25, 3 * Math.PI / 2 + 0.25, 3 * Math.PI / 2 + 0.25}};

        SchottkyDimersQuadUnitary schottkyDimers = new SchottkyDimersQuadUnitary(new SchottkyData(schottkyParamsCol), angles);
        return schottkyDimers;
    }

}
