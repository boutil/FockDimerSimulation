package dimerSim;

import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import de.jtem.riemann.theta.Theta;

public class Z2LatticeFock extends Z2Lattice{
    

    public SchottkyDimersQuad schottkyDimers;
    public ComplexVector[][] discreteAbelMap;
    // public int N, M;
    public ComplexVector Z;
    public Theta theta;
    public double[][] faceWeightsBeforeEFactors;

    // For now we assume only 4 angles: alpha-, beta-, alpha+, beta+
    public Complex[] angles;


    public Z2LatticeFock(SchottkyDimersQuad dimers, int N, int M) {
        super(N, M, 1);
        schottkyDimers = dimers;
        theta = new Theta(schottkyDimers.getPeriodMatrix());
        angles = schottkyDimers.getAngles();
        // TODO initialize this in a sensible way? 0 for now.
        Z = new ComplexVector(dimers.getNumGenerators());
        discreteAbelMap = new ComplexVector[N+2][M+2];
        flipFaceWeights = new double[N][M];
        faceWeightsBeforeEFactors = new double[N][M];
        computeDiscreteAbelMap();
        computeFaceWeights();
    }

    
    private Complex getAngle(boolean isAlpha, int k) {
        // returns the kth alpha or beta angle. Assumes total of 4 angles for now
        // Alpha sequence: [a_0^-, a_0^+, a_1^-, ...]
        // Beta sequence:  [b_0^+, b_1^-, b_1^+, ...]
        if (isAlpha) {
            // Here we use floormod instead of % to ensure that the number is positive.
            return angles[2 * Math.floorMod(k, 2)];
        } else {
            return angles[3 - 2 * Math.floorMod(k, 2)];
        }
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

     private void computeDiscreteAbelMap() {
         // dynamically compute the discrete Abel maps
         for(int i = 0; i < discreteAbelMap.length; i++) {
             if(i == 0){
                 discreteAbelMap[0][0] = Z;
             }
             else{
                 discreteAbelMap[i][0] = discreteAbelMap[i - 1][0].plus(discreteMapUpdateStep(i, 0, 0));
             }
             for (int j = 1; j < discreteAbelMap.length; j++) {
                 discreteAbelMap[i][j] = discreteAbelMap[i][j-1].plus(discreteMapUpdateStep(i, j, 1));
             }
         }
     }

    private ComplexVector discreteMapUpdateStep(int i, int j, int direction) {
        // eta(i, j) - eta(i - 1, j) if direction == 0
        // eta(i, j) - eta(i, j - 1) if direction == 1
        // Handles the angles
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
            return v1.minus(v2).divide(new Complex(0, 2*Math.PI));
        } else {
            return v2.minus(v1).divide(new Complex(0, 2*Math.PI));
        }
    }

    public void computeFaceWeights() {
        Complex[] thetaSums = new Complex[4];
        Complex[] factors = new Complex[4];
        for (int i = 0; i < factors.length; i++) {
            thetaSums[i] = new Complex();
            factors[i] = new Complex();
        }
        // crossRatio of theta of the surrounding faces (and the E differentials)
        for (int i = 1; i < N + 1; i++) {
            for (int j = 1; j < M + 1; j++) {
                theta.theta(discreteAbelMap[i-1][j], factors[0], thetaSums[0]);
                theta.theta(discreteAbelMap[i+1][j], factors[2], thetaSums[2]);
                theta.theta(discreteAbelMap[i][j-1], factors[1], thetaSums[1]);
                theta.theta(discreteAbelMap[i][j+1], factors[3], thetaSums[3]);
                // stable way of computing the crossratio treats exp part and sum part separately like this.
                Complex factorsCrossSum = factors[0].plus(factors[2]).minus(factors[1]).minus(factors[3]);
                Complex crossRatio = Complex.exp(factorsCrossSum).times(thetaSums[0]).times(thetaSums[2]).divide(thetaSums[1]).divide(thetaSums[3]);
                // Complex crossRatio = theta.theta(discreteAbelMap[i+1][j]).times(theta.theta(discreteAbelMap[i-1][j])).divide(
                //     theta.theta(discreteAbelMap[i][j+1]).times(theta.theta(discreteAbelMap[i][j-1]))
                // );
                Complex[] faceAngles = getAnglesOfFace(i, j);
                Complex quotiontOfEs1 = schottkyDimers.abelianIntegralOf3rdKind(faceAngles[1], faceAngles[0], faceAngles[2]).exp();
                Complex quotiontOfEs2 = schottkyDimers.abelianIntegralOf3rdKind(faceAngles[3], faceAngles[2], faceAngles[0]).exp();
                Complex Equotients = quotiontOfEs1.times(quotiontOfEs2);
                crossRatio = crossRatio.times(quotiontOfEs1).times(quotiontOfEs2);

                flipFaceWeights[i-1][j-1] = -crossRatio.re; // Should be like this. Check that these are indeed real first!  
                flipFaceWeights[i-1][j-1] *= ((i+j)%2 == 0) ? volumeConstraint : 1/volumeConstraint;
            }
        }
    }

}
