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

    public Z2Lattice(int N, int M) {
        flipFaceWeights = new double[N][M];
        this.N = N;
        this.M = M;
        for(double[] row : flipFaceWeights){
            Arrays.fill(row, 1);
        }
    }

    public double getUnflippedFaceWeight(int i, int j) {
        if ((i + j)%2 == 0) {
            return flipFaceWeights[i][j];
        } else {
            return 1/flipFaceWeights[i][j];
        }
    }
}
