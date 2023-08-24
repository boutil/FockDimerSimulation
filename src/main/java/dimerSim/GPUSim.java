package dimerSim;

import java.util.Random;

import jcuda.Pointer;
import jcuda.runtime.JCuda;
import lattices.Lattice;

public class GPUSim {
    // class that will run simulation steps on CUDA compatible GPU.

    Lattice lattice;

    byte[][] faceStates;
    boolean[][] insideBoundary;
    
    Random rand;


    public GPUSim(Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        this.lattice = lattice;
        this.faceStates = faceStates;
        this.insideBoundary = insideBoundary;
        rand = new Random();
    }

    private void initCuda() {
        // Copy memory and load kernel functions
        int N = lattice.N;

        // separate faceWeights into two arrays.
        double[] faceWeightsEven = new double[N * (N/2)];
        double[] faceWeightsOdd = new double[N * (N/2)];

        Pointer faceWeightsEvenPointer = new Pointer().to(faceWeightsEven);

        JCuda.cudaMemcpy2D(null, N, null, N, N, N, N);

    }



    public void simulate(int steps) {
        initCuda();

        for (int i = 0; i < steps; i++) {
            int parity = rand.nextInt(2);
            if (parity == 0) {

            } else {

            }
        }
    }

}
