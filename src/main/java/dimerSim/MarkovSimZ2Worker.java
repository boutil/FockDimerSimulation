package dimerSim;

import java.util.Random;

public class MarkovSimZ2Worker extends Thread{

    private boolean running;
    public volatile boolean doStep;
    private volatile int stepParity;
    public volatile boolean doCleanup;
    private volatile int cleanupParity;

    private MarkovSimZ2 sim;
    // The rows this worker is responsible for.
    private int[] rowIndices;
    Random rand;
    private byte ns = 0b1010;

    public MarkovSimZ2Worker(MarkovSimZ2 sim, int[] rowIndices) {
        this.sim = sim;
        this.rowIndices = rowIndices;
        rand = new Random();
    }

    @Override
    public void run() {
        running = true;
        while(running) {
            if(doStep) {
                markovStep();
                // System.out.println("done with step in thread " + rowIndices[0] / 70);
                doStep = false;
            }
            if(doCleanup) {
                cleanupStep();
                // System.out.println("cleanup done in thread  " + rowIndices[0] / 70);
                doCleanup = false;
            }
            if (Thread.interrupted()) {
                return;
            }
        }
    }

    public void requestStep(int parity) {
        this.stepParity = parity;
        this.doStep = true;
    }
    
    public void requestCleanup(int parity) {
        this.cleanupParity = parity;
        this.doCleanup = true;
    }

    public void stopRunning() {
        running = false;
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
