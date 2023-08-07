package dimerSim;

import java.util.Random;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

public class MarkovSimWorker extends Thread{

    protected volatile boolean running;
    public volatile boolean doStep;
    protected volatile int stepParity;
    public volatile boolean doCleanup;
    protected volatile int cleanupParity;

    protected MarkovSim sim;
    // The rows this worker is responsible for.
    protected int[] rowIndices;
    Random rand;

    public MarkovSimWorker(MarkovSim sim, int[] rowIndices) {
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
        // requests a step on all faces of given parity.
        this.stepParity = parity;
        this.doStep = true;
    }
    
    public void requestCleanup(int parity) {
        // Should request a cleanup on all faces that are of different parity than parity.
        this.cleanupParity = parity;
        this.doCleanup = true;
    }

    public void stopRunning() {
        running = false;
    }

    public void cleanupStep() {
        throw new NotImplementedException();
    }

    public void markovStep() {
        throw new NotImplementedException();
    }
}
