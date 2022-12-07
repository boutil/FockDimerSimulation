package dimerSim;

import org.checkerframework.checker.units.qual.K;

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
    public Complex[][] faceWeights;

    // For now we assume only 4 angles: alpha-, beta-, alpha+, beta+
    public Complex[] angles;

    public Z2LatticeFock(SchottkyDimersQuad dimers, int N, int M) {
        super(N, M);
        schottkyDimers = dimers;
        theta = new Theta(schottkyDimers.getPeriodMatrix());
        angles = schottkyDimers.getAngles();
        // TODO initialize this in a sensible way? 0 for now.
        Z = new ComplexVector(dimers.getNumGenerators());
        discreteAbelMap = new ComplexVector[N+2][M+2];
        flipFaceWeights = new double[N][M];
        faceWeights = new Complex[N][M];
        computeDiscreteAbelMap();
        computeFaceWeights();
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

    private Complex getAngle(boolean isAlpha, int k) {
        // returns the kth alpha or beta angle. Assumes total of 4 angles for now
        if (isAlpha) {
            // Here we use floormod instead of % to ensure that the number is positive.
            return angles[2 * Math.floorMod(k, 2)];
        } else {
            return angles[3 - 2 * Math.floorMod(k, 2)];
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
        if (direction == 1){
            alpha = getAngle(true, i + j);
            beta = getAngle(false, j - i);
        } else {
            alpha = getAngle(true, i + j + 1);
            beta = getAngle(false, j - i);
        }
        schottkyDimers.abelMap(v1, alpha);
        schottkyDimers.abelMap(v2, beta);
        if((i+j)%2 == 0) {
            return v1.minus(v2);
        } else {
            return v2.minus(v1);
        }
    }

    // private Complex[] getAnglesAroundFace(int i, int j) {

    // }

    public void computeFaceWeights() {
        // crossRatio of theta of the surrounding faces (and the E differentials)
        for (int i = 1; i < N + 1; i++) {
            for (int j = 1; j < M + 1; j++) {
                // faceWeights[i - 1][j - 1] = crossRatio
                Complex crossRatio = theta.theta(discreteAbelMap[i+1][j]).times(theta.theta(discreteAbelMap[i-1][j])).divide(
                    theta.theta(discreteAbelMap[i][j+1]).times(theta.theta(discreteAbelMap[i][j-1]))
                );
                faceWeights[i-1][j-1] = crossRatio;
                // TODO include the E factor crossratio here.
                flipFaceWeights[i-1][j-1] = crossRatio.re; // Should be like this. Check that these are indeed real first!  
            }
        }
    }

}
