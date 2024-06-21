package lattices;

import de.jtem.blas.ComplexVector;
import de.jtem.blas.IntegerVector;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import de.jtem.riemann.theta.Theta;
import de.jtem.riemann.theta.ThetaWithChar;

public class Z2LatticeFock extends Z2Lattice{
    

    public transient SchottkyDimers schottkyDimers;
    public ComplexVector[][] discreteAbelMap;
    // public int N, M;
    public ComplexVector Z;
    public transient Theta theta;

    private transient ThetaWithChar thetaWithChar;

    // angles of the form [[alpha^-_0, alpha^-_1, ...][beta^-_0, beta^-_1, ...][alpha^+_0, alpha^+_1, ...][beta^+_0, beta^+_1, ...]]
    public Complex[][] angles;
    public ComplexVector[][] faceAngleAbelMaps = new ComplexVector[4][];

    // We impose that the number of all angle types must be equal for simplicity of periodicity concerns.
    private int numAlphas;
    public ComplexVector[][] abelIncrementsRight;
    public ComplexVector[][] abelIncrementsTop;


    public Z2LatticeFock(SchottkyDimers dimers, int N, int M) {
        super(N, M, 1);
        schottkyDimers = dimers;
        // theta = new Theta(schottkyDimers.getPeriodMatrix());
        theta = new Theta(schottkyDimers.getPeriodMatrix(), 1e-7, false);
        // Create a theta function with odd characteristic (1, 0, 0, ...).
        thetaWithChar = new ThetaWithChar(schottkyDimers.getPeriodMatrix(), 1e-7, false);
        int[] oddAlpha = new int[schottkyDimers.getNumGenerators()];
        oddAlpha[0] = 1;
        IntegerVector e1 = new IntegerVector(oddAlpha);
        thetaWithChar.setAlpha(e1);
        thetaWithChar.setBeta(e1);

        angles = schottkyDimers.getAngles();
        // we assume same amount of all angle types for now.
        numAlphas = angles[0].length + angles[2].length;
        abelIncrementsRight = new ComplexVector[numAlphas][numAlphas];
        abelIncrementsTop = new ComplexVector[numAlphas][numAlphas];
        for (int i = 0; i < angles.length; i++) {
            faceAngleAbelMaps[i] = new ComplexVector[angles[i].length];
            for (int k = 0; k < angles[i].length; k++) {
                faceAngleAbelMaps[i][k] = new ComplexVector(schottkyDimers.getNumGenerators());
                schottkyDimers.abelMap(faceAngleAbelMaps[i][k], angles[i][k]);
            }
        }
        // TODO initialize this in a sensible way? 0 for now.
        Z = new ComplexVector(dimers.getNumGenerators(), 0, 0);
        discreteAbelMap = new ComplexVector[N+2][M+2];
        flipFaceWeights = new double[N][M];

        for (int i = 0; i < factors.length; i++) {
            thetaSums[i] = new Complex();
            factors[i] = new Complex();
        }

        computeDiscreteAbelMap();
        computeFaceWeights();
        // checkDiscreteAbelMapCorrectness();
    }
    
    private Complex getAngle(boolean isAlpha, int k) {
        // Alpha sequence: [a_0^-, a_0^+, a_1^-, ...]
        // Beta sequence:  [b_0^+, b_1^-, b_1^+, ...]
        int index = isAlpha ? 2 * Math.floorMod(k, 2) : 3 - 2 * Math.floorMod(k, 2);
        int length = angles[index].length;
        return angles[index][Math.floorMod(k, 2 * length) / 2];
    }

    private ComplexVector getAngleAbelMap(boolean isAlpha, int k) {
        // Alpha sequence: [a_0^-, a_0^+, a_1^-, ...]
        // Beta sequence:  [b_0^+, b_1^-, b_1^+, ...]
        int index = isAlpha ? 2 * Math.floorMod(k, 2) : 3 - 2 * Math.floorMod(k, 2);
        int length = angles[index].length;
        return faceAngleAbelMaps[index][Math.floorMod(k, 2 * length) / 2];
    }

    private Complex[] getAnglesOfFace(int i, int j) {
        // returns list of  angles [alphaNW, betaSW, alphaSE, betaNE] that are associated to the given quad.
        Complex[] faceAngles = new Complex[4];
        faceAngles[0] = getAngle(true, i + j);
        faceAngles[1] = getAngle(false, j - i);
        faceAngles[2] = getAngle(true, i + j + 1);
        faceAngles[3] = getAngle(false, j - i + 1);
        return faceAngles;
     }

    private ComplexVector[] getAngleAbelMapsOfFace(int i, int j) {
    // returns list of  angle abel maps [alphaNW, betaSW, alphaSE, betaNE] that are associated to the given quad.
    ComplexVector[] faceAngles = new ComplexVector[4];
    faceAngles[0] = getAngleAbelMap(true, i + j);
    faceAngles[1] = getAngleAbelMap(false, j - i);
    faceAngles[2] = getAngleAbelMap(true, i + j + 1);
    faceAngles[3] = getAngleAbelMap(false , j - i + 1);
    return faceAngles;
    }

    private void computeDiscreteAbelMap() {
        // dynamically compute the discrete Abel maps
        computeMapUpdateSteps();
        for(int i = 0; i < discreteAbelMap.length; i++) {
            if(i == 0){
                discreteAbelMap[0][0] = Z;
            }
            else{
                discreteAbelMap[i][0] = discreteAbelMap[i - 1][0].plus(abelIncrementsRight[i % numAlphas][0]);
            }
            for (int j = 1; j < discreteAbelMap.length; j++) {
                discreteAbelMap[i][j] = discreteAbelMap[i][j-1].plus(abelIncrementsTop[i % numAlphas][j % numAlphas]);
            }
        }
    }

    private ComplexVector computeDiscreteUpdateStep(int i, int j, int direction) {
        ComplexVector v1 = new ComplexVector(schottkyDimers.getNumGenerators());
        ComplexVector v2 = new ComplexVector(schottkyDimers.getNumGenerators());
        Complex alpha;
        Complex beta;
        Complex[] faceAngles = getAnglesOfFace(i, j);
        if (direction == 1){
            alpha = faceAngles[0];
            beta = faceAngles[1];
        } else {
            alpha = faceAngles[2];
            beta = faceAngles[1];
        }
        schottkyDimers.abelMap(v1, alpha);
        schottkyDimers.abelMap(v2, beta);
        if((i+j)%2 == 0) {
            return v1.minus(v2);
        } else {
            return v2.minus(v1);
        }
    }
    // This needs to be bigger for not only 4 angles case.
    private void computeMapUpdateSteps() {
        // eta(i, j) - eta(i - 1, j) if direction == 0
        // eta(i, j) - eta(i, j - 1) if direction == 1
        // Handles the angles
        for (int i = 0; i < numAlphas; i++) {
            for (int j = 0; j < numAlphas; j++) {
                abelIncrementsRight[i][j] = computeDiscreteUpdateStep(i, j, 0);
                abelIncrementsTop[i][j] = computeDiscreteUpdateStep(i, j, 1);
            }
        }
    }
    
    Complex[] thetaSums = new Complex[4];
    Complex[] factors = new Complex[4];

    public void computeFaceWeights() {
        // crossRatio of theta of the surrounding faces (and the E differentials)
        for (int i = 1; i < N + 1; i++) {
            for (int j = 1; j < M + 1; j++) {
                theta.theta(discreteAbelMap[i-1][j], factors[0], thetaSums[0]);
                theta.theta(discreteAbelMap[i+1][j], factors[2], thetaSums[2]);
                theta.theta(discreteAbelMap[i][j-1], factors[1], thetaSums[1]);
                theta.theta(discreteAbelMap[i][j+1], factors[3], thetaSums[3]);
                // numerically stable way of computing the crossratio treats exp part and sum part separately like this.
                Complex factorsCrossSum = factors[0].plus(factors[2]).minus(factors[1]).minus(factors[3]);
                Complex crossRatio = Complex.exp(factorsCrossSum).times(thetaSums[0]).times(thetaSums[2]).divide(thetaSums[1]).divide(thetaSums[3]);
                // Complex crossRatio = theta.theta(discreteAbelMap[i+1][j]).times(theta.theta(discreteAbelMap[i-1][j])).divide(
                //     theta.theta(discreteAbelMap[i][j+1]).times(theta.theta(discreteAbelMap[i][j-1]))
                // );
                
                crossRatio.assignTimes(getECrossratio(i, j));

                double maxImagPart = 0.001;
                if(Math.abs(crossRatio.im) > maxImagPart ) {
                    System.out.println("crossRatio imaginary part too big: " + crossRatio);
                }
                if(crossRatio.re > 0 ) {
                    System.out.println("weight has wrong sign.");
                }

                flipFaceWeights[i-1][j-1] = -crossRatio.re; // Should be like this. Check that these are indeed real first!  
                flipFaceWeights[i-1][j-1] *= ((i+j)%2 == 0) ? faceWeightMultiplier : 1/faceWeightMultiplier;
            }
        }
    }

    public Complex getThetaCrossRatio(ComplexVector Z, ComplexVector topStep, ComplexVector rightStep) {
        // Evaluates the theta of theta(Z - topStep/2) * theta(Z + topStep/2) / (theta(Z - rightStep/2) * theta(Z + rightStep/2)) 
        // For visualization and sanity check purposes.
        theta.theta(Z.minus(rightStep.divide(2)), factors[0], thetaSums[0]);
        theta.theta(Z.plus(rightStep.divide(2)), factors[2], thetaSums[2]);
        theta.theta(Z.minus(topStep.divide(2)), factors[1], thetaSums[1]);
        theta.theta(Z.plus(topStep.divide(2)), factors[3], thetaSums[3]);
        // numerically stable way of computing the crossratio treats exp part and sum part separately like this.
        Complex factorsCrossSum = factors[0].plus(factors[2]).minus(factors[1]).minus(factors[3]);
        Complex crossRatio = Complex.exp(factorsCrossSum).times(thetaSums[0]).times(thetaSums[2]).divide(thetaSums[1]).divide(thetaSums[3]);
        return crossRatio;
    }

    private Complex getECrossratio(int i, int j) {
        // Complex[] faceAngles = getAnglesOfFace(i, j);

        // // Method 1 of computing this:
        // Complex quotiontOfEs1 = schottkyDimers.abelianIntegralOf3rdKind(faceAngles[1], faceAngles[0], faceAngles[2]).exp();
        // Complex quotiontOfEs2 = schottkyDimers.abelianIntegralOf3rdKind(faceAngles[3], faceAngles[2], faceAngles[0]).exp();
        // Complex eCrossRatio = quotiontOfEs1.times(quotiontOfEs2);
        // Method 2: through odd theta functions:
        // TODO: switch to numerically stable version here too!
        ComplexVector[] faceAngleAbelMaps = getAngleAbelMapsOfFace(i, j);
        Complex thetaCrossRatio = thetaWithChar.theta(faceAngleAbelMaps[2].minus(faceAngleAbelMaps[1]));
        thetaCrossRatio.assignTimes(thetaWithChar.theta(faceAngleAbelMaps[0].minus(faceAngleAbelMaps[3])));
        thetaCrossRatio.assignDivide(thetaWithChar.theta(faceAngleAbelMaps[1].minus(faceAngleAbelMaps[0])));
        thetaCrossRatio.assignDivide(thetaWithChar.theta(faceAngleAbelMaps[3].minus(faceAngleAbelMaps[2])));

        // // Compare the two results for sanity:
        // if (! thetaCrossRatio.equals(eCrossRatio, 1e-5)) {
        //     System.out.println("E CrossRatio computation methods do not agree!");
        // }

        return thetaCrossRatio.invert();
    }

}
