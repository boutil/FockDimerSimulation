package dimerSim;

import java.io.Serializable;
import java.util.List;
import java.util.Random;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

import lattices.Lattice;

public class MarkovSim implements Serializable{
    public Lattice lattice;
    // states of each face. Encodes a dimer configuration. Can compute a height function from this.
    public byte[][] faceStates;

    public int currentVolume;

    public int[][] heightFunction;
    // Indicates for each face whether it should be updated.
    public boolean[][] insideBoundary;
    
    protected Random rand;

    protected MarkovSimWorker[] markovWorkers;

    public int maxParity;

    public int numThreads = 8;

    List<Index> upFlippableIndices;
    List<Index> downFlippableIndices;
    
    // Only makes sense to not set to 1 in case of uniform 1 face weights.
    public double acceptanceRatioConstant = 0.9999;

    public MarkovSim(Lattice lattice) {
        this.lattice = lattice;
        // initialize height function in a sensible way
        faceStates = new byte[lattice.N][lattice.M];
        insideBoundary = new boolean[lattice.N][lattice.M];
    }

    public MarkovSim(Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        this.lattice = lattice;
        this.faceStates = faceStates;
        this.insideBoundary = insideBoundary;
    }

    protected void init() {
        throw new NotImplementedException();
    }

    public void setLattice(Lattice newLattice) {
        // Used to change weights for a given lattice.
        if (lattice.N != newLattice.N || lattice.M != newLattice.M) {
            System.out.println("cannot swap lattices. Dimensions do not match.");
            return;
        }
        lattice = newLattice;
    }

    protected void createWorkers(int numThreads) {

    }

    protected void markovStep(int parity) {
        for (MarkovSimWorker worker : markovWorkers) {
            worker.requestStep(parity);
        }
        while(true) {
            boolean working = false;
            for (MarkovSimWorker worker : markovWorkers) {
                working |= worker.doStep;
            }
            if(!working) {
                break;
            }
        }
        for(int paritySummand = 1; paritySummand < maxParity; paritySummand++){
            for (MarkovSimWorker worker : markovWorkers) {
                worker.requestCleanup((parity + paritySummand) % maxParity);
            }
            while(true) {
                boolean working = false;
                for (MarkovSimWorker worker : markovWorkers) {
                    working |= worker.doCleanup;
                }
                if(!working) {
                    break;
                }
            }
        }
    }

    public boolean isDimer(Index coords, int direction) {
        // returns true if there is a dimer in the given direction from the coords face.
        return ((faceStates[coords.x][coords.y] >> direction) & 1) != 0;
    }

    public void flipFace(Index ind) {
        throw new NotImplementedException();
    }

    public void flipFaceExclusive(int i, int j) {
        throw new NotImplementedException();
    }

    public void consolidateFaceStateFromNeighbors(int i, int j) {
        throw new NotImplementedException();
    }

    public void flipFace(int i, int j) {
        flipFace(new Index(i,j));
    }

    public int flippableDirection(Index ind) {
        return flippableDirection(ind.x, ind.y);
    }

    public int flippableDirection(int i, int j) {
        // returns 1 if flippable in + direction, -1 if in - direction and 0 if not flippable
        throw new NotImplementedException();
    }

    protected int computeVolume() {
        // computes volume enclosed by the height function. That's just the integral
        int vol = 0;
        for (int i = 0; i < faceStates.length; i++) {
            for (int j = 0; j < faceStates.length; j++) {
                vol += heightFunction[i][j];
            }
        }
        return vol;
    }

    public void simulate(int numSteps) {
        createWorkers(numThreads);
        for (MarkovSimWorker worker : markovWorkers) {
            worker.start();
        }
        int reportFreq = numSteps/10;
        long time = System.currentTimeMillis();
        for (int i = 0; i < numSteps; i++) {
            if ((i % reportFreq == 0)) {
                double timeForAvg1000Steps = Math.max(((double)(System.currentTimeMillis() - time)) * 1000 / reportFreq, 1000);
                System.out.println("Done with " + i + " steps." + " Average time per 1000 markovSteps: " + (int) timeForAvg1000Steps + ". Time left: " + (int)((numSteps - i) * timeForAvg1000Steps / 1000000) + " seconds.");
                time = System.currentTimeMillis();
            }
            // Choose parity at random at each step. Can probably just alternate too? 
            // Do we need to take weights into account here?
            markovStep(rand.nextInt(maxParity));
        }
        for (MarkovSimWorker worker : markovWorkers) {
            worker.stopRunning();
        }
    }

    public Integer getHeight(int i, int j) {
        if (insideBoundary[i][j]) {
            return heightFunction[i][j];
        } else{
            return 0;
        }
    }


    protected void computeHeightFunctionFromFaceStates() {
        // computes the height function associated to the current face states.
        // Note that this only produces an actual height function if the face states encode a perfect matching! - (Could check for this - maybe todo)
        throw new NotImplementedException();
    }
}
