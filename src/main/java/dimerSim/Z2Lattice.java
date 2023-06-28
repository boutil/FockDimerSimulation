package dimerSim;

import java.io.Serializable;
import java.util.Arrays;

import de.jtem.mfc.field.Complex;

public class Z2Lattice implements Serializable{
    // The flipFaceWeights here are defined as the crossratio of (N * S)/(E * W). 
    // !!CAREFUL!! Note that this is is not the same thing as the usual faceWeight definition. It is the same for white faces but the inverse for black faces!
    //  However it allows to do the flip easily in the simulation without considering many cases.
    public double[][] flipFaceWeights;
    public int N, M;

    public double faceWeightMultiplier;

    public Z2Lattice(int N, int M, double faceWeightMultiplier) {
        flipFaceWeights = new double[N][M];
        this.N = N;
        this.M = M;

        this.faceWeightMultiplier = faceWeightMultiplier;
        for(double[] row : flipFaceWeights){
            Arrays.fill(row, 1);
        }
        for (int i = 0; i < flipFaceWeights.length; i++) {
            for (int j = 0; j < flipFaceWeights[i].length; j++) {
                flipFaceWeights[i][j] = ((i+j)%2 == 0) ? faceWeightMultiplier : 1/faceWeightMultiplier;
                // flipFaceWeights[i][j] = volumeConstraint;
            }
        }
    }

    public Z2Lattice(int N, int M) {
        this(N, M, 1);
    }

    public double getUnflippedFaceWeight(int i, int j) {
        if ((i + j)%2 == 0) {
            return flipFaceWeights[i][j];
        } else {
            return 1/flipFaceWeights[i][j];
        }
    }
}
