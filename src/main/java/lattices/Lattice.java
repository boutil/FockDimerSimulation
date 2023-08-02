package lattices;

import java.io.Serializable;

public abstract class Lattice implements Serializable{
    public double[][] flipFaceWeights;
    public int N, M;
    public double faceWeightMultiplier;

    public Lattice(int N, int M, double faceWeightMultiplier) {
        flipFaceWeights = new double[N][M];
        this.N = N;
        this.M = M;

        this.faceWeightMultiplier = faceWeightMultiplier;
    }

    public Lattice(int N, int M) {
        this(N, M, 1);
    }

    public abstract double getUnflippedFaceWeight(int i, int j);
}
