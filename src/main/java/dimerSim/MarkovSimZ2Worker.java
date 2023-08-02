package dimerSim;

import java.util.Random;

public class MarkovSimZ2Worker extends MarkovSimWorker{

    private byte ns = 0b1010;

    public MarkovSimZ2Worker(MarkovSimZ2 sim, int[] rowIndices) {
        super(sim, rowIndices);
    }

    public void cleanupStep() {
        int parity = cleanupParity;
        for (int index = 0; index < rowIndices.length; index++) {
            int i = rowIndices[index];
            for (int j = 0; j < sim.faceStates[i].length; j++) {
                if (((i + j) % 2 == parity)) {
                    sim.consolidateFaceStateFromNeighbors(i, j);
                }
            }
        }
    }

    public void markovStep() {
        int parity = stepParity;
        for (int index = 0; index < rowIndices.length; index++) {
            int i = rowIndices[index];
            for (int j = 0; j < sim.faceStates[i].length; j++) {
                if (((i + j) % 2 == parity) && sim.insideBoundary[i][j]) {
                    double upProp = sim.lattice.flipFaceWeights[i][j];
                    double prop = (sim.faceStates[i][j] == ns) ? upProp : 1/upProp;
                    prop *= sim.acceptanceRatioConstant;
                    if (rand.nextDouble() < prop) {
                        sim.flipFaceExclusive(i, j);
                    }
                }
            }
        }
    }
    
}
