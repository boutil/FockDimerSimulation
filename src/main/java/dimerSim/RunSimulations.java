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
        // int numSteps = (int) 1e2;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Hexagon/03_convergenceToUniform/";
        // int size = 600;
        // int numSteps = (int)1e6;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Hexagon/04_growingStarfishHoles/mu06/";
        // int size = 600;
        // int numSteps = (int)1e3;
        // int interval = 0;
        // int numTimes = 1;

        String experimentName = "Hexagon/05_g0anglesOn2Sides/mu00000001/";
        int size = 600;
        int numSteps = (int)1e2;
        int interval = 0;
        int numTimes = 1;

        // --------------------------------------------------

        // String experimentName = "Aztec/01_LargeHoleG1/mu03/";
        // int size = 600;
        // int numSteps = (int) 5e6;
        // int interval = (int) 0;
        // int numTimes = 4;
        
        // String experimentName = "Aztec/02_LargeHole2Angles/mu05/";
        // int size = 600;
        // int numSteps = (int)1e3;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/03_G1Unitary3Angles/mu08/";
        // int size = 601;
        // int numSteps = (int) 1e3;
        // int interval = 0;
        // int numTimes = 1;

        // String experimentName = "Aztec/04_G2LargeHoles3Angles/mu00003/";
        // int size = 601;
        // int numSteps = (int) 1e3;
        // int interval = 0;
        // int numTimes = 1;

        // String experimentName = "Aztec/05_growing1LargeHoleUnitary/mu1/";
        // int size = 600;
        // int numSteps = (int) 1e3;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/06_growing2LargeHolesFromUnitary/mu01/";
        // int size = 601;
        // int numSteps = (int) 1e3;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/07_growing2LargeHolesFromUnitarySymmetric/mu03asd/";
        // int size = 600;
        // int numSteps = (int) 10;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/08_growing2LargeHolesFromUnitarySymmetric1Angle/mu03/";
        // int size = 600;
        // int numSteps = (int) 1e3;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/09_growing2LargeHolesFromUnitaryWithAngles/mu01/";
        // int size = 600;
        // int numSteps = (int) 1e3;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/10_growing2LargeHolesFromUnitarySymmetric1Angle/mu03/";
        // int size = 600;
        // int numSteps = (int) 1e7;
        // int interval = 0;
        // int numTimes = 5;

        // String experimentName = "Aztec/11_genus0_3Angles/mu000000001/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 1;

        // String experimentName = "Aztec/12_genus0_3Angles2Sides/mu000000001/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 1;

        // String experimentName = "Aztec/13_genus0_3AnglesAllSides/mu000000001/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 5;
        
        // String experimentName = "Aztec/14_genus0_2AnglesAllSides/mu000000001/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 3;

        // String experimentName = "Aztec/15_genus0_3Angles2Sides/mu000000001/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 5;

        // String experimentName = "Aztec/16_genus0_3StrongEllipse/mu000000001/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 5;
        
        // String experimentName = "Aztec/17_genus0_4Anglesasd/mu000000001/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/18_genus1_central/mu2/";
        // int size = 600;
        // int numSteps = (int) 1e6;
        // int interval = 0;
        // int numTimes = 1;
        
        // String experimentName = "Aztec/20_genus1_3Angles/mu0000001asd/";
        // int size = 600;
        // int numSteps = (int) 1e2;
        // int interval = 0;
        // int numTimes = 1;

        System.out.println("Running experiment " + experimentName + " of size " + size);
        
        String folderName = baseFolder + experimentName;
        
        // SchottkyDimers schottkyDimers = createSchottkyQuad();
        // // SchottkyDimers schottkyDimers = createSchottkyHex();

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
        // // Starfish:
        // Complex A = new Complex(Math.cos(Math.PI / 3 + 0.01), Math.sin(Math.PI / 3 + 0.01)).times(0.5);
        // double d = 0.07;
        // double[][] angles = {{d, Math.PI / 3 - d}, {Math.PI / 3 + d, 2 * Math.PI / 3 - d}, {2 * Math.PI / 3 + d, 3 * Math.PI / 3 - d}};

        // One bubble with angles
        Complex A = new Complex(Math.cos(2 * Math.PI / 3 + 0.1), Math.sin(2 * Math.PI / 3 + 0.1)).times(0.5);
        // double[][] angles = {{0, 0}, {Math.PI / 3 - 0.6, Math.PI / 3 + 0.6}, {2 * Math.PI / 3, 2 * Math.PI / 3}};
        
        Complex B = A.invert().conjugate();
        // double[] schottkyParams = {A.re, A.im, B.re, B.im, 0.05, 0};
        
        // 05_g0AnglesOn2Sides
        double[][] angles = {{0, 0}, {Math.PI / 3, Math.PI / 3}, {2 * Math.PI / 3 - 0.9, 2 * Math.PI / 3 + 0.9}};
        double[] schottkyParams = {A.re, A.im, B.re, B.im, 0.0000000001, 0};

        // double[] sizes = {1, 1, 1};
        double[] sizes = {1, 0.8, 0.6};

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
        
        // // 02_G1LargeHole2Angles
        // double theta1 = -Math.PI/4;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.2);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.00000001, 0};
        // double[][] angles = {{0, 0}, {Math.PI / 2 - 0.9, Math.PI / 2 + 0.4}, {Math.PI, Math.PI}, {3 * Math.PI / 2, 3 * Math.PI / 2}};

        // // 04_G2LargeHoles3Angles
        // double[] schottkyParams = {-0.54, 0.69, -0.54, -0.69, 0.00003, 0, 0.69, 0.89, 0.69, -0.89, 0.00003, 0};
        // // double[] schottkyParams = {-2.4, 0.69, -2.4, -0.69, 0.03, 0, 0.69, 0.89, 0.69, -0.89, 0.03, 0};
        // double[][] angles = {{-2.4, -2.4, -2.4}, {-1.2, -0.4, 0.1}, {0.4, 0.4, 0.4}, {0.8, 2.4, 15}};

        // // // 04_G2LargeHoles3Angles
        // double[] schottkyParams = {-0.5, 0.7, -0.5, -0.7, 0.00000000001, 0};
        // double[][] angles = {{-2.4, -2.4, -2.4}, {-1.2, -0.4, 0.1}, {0.15, 0.4, 0.4}, {0.6, 2, 15}};
        
        // // 16_genus0_3StrongEllipse
        // double[] schottkyParams = {0, 1, 0, -1, 0.00000000001, 0};
        // double[][] angles = {{-2.4}, {-0.4}, {1.5}, {2.4}};

        // 17_genus0_4Angles
        // double[] schottkyParams = {0, 1, 0, -1, 0.00000000001, 0};
        // double[][] angles = {{-2.4, -2.4, -2.4, -2.4}, {-1.2, -0.8, -0.45, -0.05}, {0.05, 0.45, 0.8, 1.2}, {2.4, 2.4, 2.4, 2.4}};

        // SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(new SchottkyData(schottkyParams), angles);
        // return schottkyDimers;

        // // 03_G1Unitary3Angles
        // // double theta1 = Math.PI/2;
        // double theta1 = Math.PI + 0.01;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.5);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.05, 0};
        // double[][] angles = {{-1, 0, 1}, {Math.PI / 2, Math.PI / 2, Math.PI/2}, {Math.PI - 1, Math.PI, Math.PI + 1}, {3 * Math.PI / 2, 3 * Math.PI / 2, 3 * Math.PI / 2}};

        
        // // G1Unitary
        // double theta = Math.PI / 4;
        // Complex A = new Complex(Math.cos(theta), Math.sin(theta)).times(0.3);
        // Complex ARefl = A.invert().conjugate();
        // double[] schottkyParamsCol = {A.re, A.im, ARefl.re, ARefl.im, 0.03, 0};
        // double[][] angles = {{0, 0}, {Math.PI / 2 - 0.5, Math.PI / 2 + 0.5}, {Math.PI - 0.5, Math.PI + 0.5}, {3 * Math.PI / 2, 3 * Math.PI / 2}};

        // // 08_G2Unitary
        // double theta1 = Math.PI / 4;
        // // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.35);
        // Complex A1Refl = A1.invert().conjugate();
        // double theta2 = 5 * Math.PI / 4 + 0.1;
        // Complex A2 = new Complex(Math.cos(theta2), Math.sin(theta2)).times(0.35);
        // Complex A2Refl = A2.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.03, 0, A2.re, A2.im, A2Refl.re, A2Refl.im, 0.03, 0};
        // double[][] angles = {{0, 0}, {Math.PI / 2 - 0.5, Math.PI / 2 + 0.5}, {Math.PI - 0.5 + 0.1, Math.PI + 0.5 + 0.2}, {3 * Math.PI / 2, 3 * Math.PI / 2}};
        
        
        // // 09_G2Unitary
        // double theta1 = Math.PI / 2 + 0.2;
        // // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.35);
        // Complex A1Refl = A1.invert().conjugate();
        // double theta2 = 3 * Math.PI / 2 + 0.2;
        // Complex A2 = new Complex(Math.cos(theta2), Math.sin(theta2)).times(0.35);
        // Complex A2Refl = A2.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.03, 0, A2.re, A2.im, A2Refl.re, A2Refl.im, 0.03, 0};
        // double[][] angles = {{-1, 0, 1}, {Math.PI / 2, Math.PI / 2, Math.PI/2}, {Math.PI - 1, Math.PI, Math.PI + 1}, {3 * Math.PI / 2, 3 * Math.PI / 2, 3 * Math.PI / 2}};
        
        // // 10_G2Unitary
        // double theta1 = Math.PI / 2 + 0.5;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.3);
        // Complex A1Refl = A1.invert().conjugate();
        // double theta2 = -0.2;
        // Complex A2 = new Complex(Math.cos(theta2), Math.sin(theta2)).times(0.3);
        // Complex A2Refl = A2.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.03, 0, A2.re, A2.im, A2Refl.re, A2Refl.im, 0.03, 0};
        // double[][] angles = {{0, 0}, {Math.PI / 2 - 0.9, Math.PI / 2 + 0.4}, {Math.PI, Math.PI}, {3 * Math.PI / 2, 3 * Math.PI / 2}};

        // 11_genus0_3Angles
        // double theta1 = Math.PI / 2 + 0.5;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.3);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000001, 0};
        // double[][] angles = {{0, 0, 0}, {Math.PI / 2 - 1, Math.PI / 2, Math.PI / 2 + 1}, {Math.PI, Math.PI, Math.PI}, {3 * Math.PI / 2, 3 * Math.PI / 2, 3 * Math.PI / 2}};

        // // 12_genus0_3Angles2Sides
        // double theta1 = Math.PI / 2 + 0.5;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.3);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000001, 0};
        // double[][] angles = {{0, 0, 0}, {Math.PI / 2 - 1, Math.PI / 2, Math.PI / 2 + 0.6}, {Math.PI - 0.6, Math.PI, Math.PI + 1}, {3 * Math.PI / 2, 3 * Math.PI / 2, 3 * Math.PI / 2}};

        // // 13_genus0_3AnglesAllSides
        // double theta1 = Math.PI / 2 + 0.5;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.3);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000001, 0};
        // double[][] angles = {{0, 0.6, 0 + 1.2}, {Math.PI / 2, Math.PI / 2 + 0.6, Math.PI / 2 + 1.2}, {Math.PI, Math.PI + 0.6, Math.PI + 1.2}, {3 * Math.PI / 2, 3 * Math.PI / 2 + 0.6, 3 * Math.PI / 2 + 1.2}};

        // // 14_genus0_2AnglesAllSides
        // double theta1 = Math.PI / 2 + 0.5;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.3);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.000000001, 0};
        // double[][] angles = {{0, 1.35}, {Math.PI / 2, Math.PI / 2 + 1.2}, {Math.PI, Math.PI + 0.6}, {3 * Math.PI / 2 - 0.3, 3 * Math.PI / 2 + 0.3}};
        
        // // 18_genus1_central
        // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.4);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.2, 0};
        // double[][] angles = {{0}, {Math.PI / 2}, {Math.PI}, {3 * Math.PI / 2}};

        // // 19_genus1_3Angles
        // double theta1 = Math.PI / 2;
        // Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.4);
        // Complex A1Refl = A1.invert().conjugate();
        // double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.0000001, 0};
        // double[][] angles = {{0, 0, 0}, {Math.PI / 2 - 1.3, Math.PI / 2, Math.PI / 2 + 1.3}, {Math.PI, Math.PI, Math.PI}, {3 * Math.PI / 2 - 1.3, 3 * Math.PI / 2, 3 * Math.PI / 2 + 1.3}};

        // 19_genus1_3Angles
        double theta1 = Math.PI / 2;
        Complex A1 = new Complex(Math.cos(theta1), Math.sin(theta1)).times(0.4);
        Complex A1Refl = A1.invert().conjugate();
        double[] schottkyParamsCol = {A1.re, A1.im, A1Refl.re, A1Refl.im, 0.0000001, 0};
        // double[][] angles = {{0, 0, 0}, {Math.PI / 2 - 1, Math.PI / 2, Math.PI / 2 + 1}, {Math.PI - 0.3, Math.PI + 0.3, Math.PI + 0.3}, {3 * Math.PI / 2 - 1, 3 * Math.PI / 2, 3 * Math.PI / 2 + 1}};
        double[][] angles = {{0, 0, 0}, {Math.PI / 2 - 1, Math.PI / 2, Math.PI / 2 + 1}, {Math.PI - 0.3, Math.PI + 0.2, Math.PI + 0.2}, {3 * Math.PI / 2 - 1, 3 * Math.PI / 2, 3 * Math.PI / 2 + 0.85}};

        SchottkyDimersQuadUnitary schottkyDimers = new SchottkyDimersQuadUnitary(new SchottkyData(schottkyParamsCol), angles);
        return schottkyDimers;
    }


}
