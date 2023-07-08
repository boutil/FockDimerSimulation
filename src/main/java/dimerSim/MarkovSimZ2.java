package dimerSim;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;
import java.util.Arrays;

public class MarkovSimZ2 implements Serializable{
    public Z2Lattice lattice;
    // states of each face. Encodes a dimer configuration. Can compute a height function from this.
    // 2^0 * W + 2^1 * S + 2^2 * E + 2^3 * N
    public byte[][] faceStates;
    private byte ns = 0b1010;
    private byte we = 0b0101;

    public int currentVolume;

    public int[][] heightFunction;
    // Indicates for each face whether it should be updated.
    public boolean[][] insideBoundary;
    
    private Random rand;

    private MarkovSimZ2Worker[] markovWorkers;

    List<Index> upFlippableIndices;
    List<Index> downFlippableIndices;
    
    // Only makes sense to not set to 1 in case of uniform 1 face weights.
    public double acceptanceRatioConstant = 0.8;

    public MarkovSimZ2(Z2Lattice lattice, boolean flat) {
        this.lattice = lattice;
        // initialize height function in a sensible way
        faceStates = new byte[lattice.N][lattice.M];
        insideBoundary = new boolean[lattice.N][lattice.M];

        init();
        
        if (flat){
            initializeFlatSquare();
        } else {
            initializeAztecDiamond();
        }
    }

    public MarkovSimZ2(Z2Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        this.lattice = lattice;
        this.faceStates = faceStates;
        this.insideBoundary = insideBoundary;
        init();
    }

    private void init() {
        heightFunction = new int[lattice.N][lattice.M];
        long seed = 42; // for reproducability
        rand = new Random(seed);
        // for clunky parallelization purposes:
        int chunkSize = 51;
        markovWorkers = new MarkovSimZ2Worker[lattice.N / chunkSize + 1];
        for (int i = 0; i < markovWorkers.length; i++) {
            markovWorkers[i] = new MarkovSimZ2Worker(this, IntStream.range(i * chunkSize, Math.min(lattice.N, (i+1) * chunkSize)).toArray());
            markovWorkers[i].start();
        }
    }

    public void setLattice(Z2Lattice newLattice) {
        // Used to change weights for a given lattice.
        if (lattice.N != newLattice.N || lattice.M != newLattice.M) {
            System.out.println("cannot swap lattices. Dimensions do not match.");
            return;
        }
        lattice = newLattice;
    }

    private void markovStepUpVol(int parity, int volumeToAchieve) {
        // Takes a step into the direction of +volume.
        int flipSign = (volumeToAchieve - currentVolume / 4 > 0) ? 1 : -1;
        for (int i = 0; i < lattice.N; i++) {
            for (int j = 0; j < lattice.M; j++) {
                if (currentVolume/4 == volumeToAchieve/4) {
                    return;
                }
                if (!insideBoundary[i][j]) {
                    continue;
                }
                if ((i + j) % 2 == parity) {
                    if (flippableDirection(i, j) == flipSign) {
                        flipFace(i, j);
                    }
                }
            }    
        }
    }

    private void markovStep(int parity) {
        for (MarkovSimZ2Worker worker : markovWorkers) {
            worker.requestStep(parity);
        }
        while(true) {
            boolean working = false;
            for (MarkovSimZ2Worker worker : markovWorkers) {
                working |= worker.doStep;
            }
            if(!working) {
                break;
            }
        }
        for (MarkovSimZ2Worker worker : markovWorkers) {
            worker.requestCleanup((parity + 1) % 2);
        }
        while(true) {
            boolean working = false;
            for (MarkovSimZ2Worker worker : markovWorkers) {
                working |= worker.doCleanup;
            }
            if(!working) {
                break;
            }
        }
        // for (int i = 0; i < faceStates.length; i++) {
        //     for (int j = 0; j < faceStates[i].length; j++) {
        //         if (((i + j) % 2 == parity) && insideBoundary[i][j]) {
        //             // if (flippableDirection(i, j) != 0) {
        //             double upProp = lattice.flipFaceWeights[i][j];
        //             double prop = (faceStates[i][j] == ns) ? upProp : 1/upProp;
        //             prop *= acceptanceRatioConstant;
        //             if (rand.nextDouble() < prop) {
        //                 flipFace(i, j);
        //             }
        //             // }
        //         }
        //     }
        // }
    }

    private void markovStepSameVol() {
        // create pairs of up and down flppable pairs and flip them all with the right probability.
        upFlippableIndices = new ArrayList<Index>();
        downFlippableIndices = new ArrayList<Index>();
        for (int i = 0; i < faceStates.length; i++) {
            for (int j = 0; j < faceStates[i].length; j++) {
                if(flippableDirection(i, j) == 1) {
                    upFlippableIndices.add(new Index(i, j));
                } else if (flippableDirection(i, j) == -1) {
                    downFlippableIndices.add(new Index(i, j));
                }
            }
        }
        Collections.shuffle(upFlippableIndices);
        Collections.shuffle(downFlippableIndices);
        for (int i = 0; i < Math.min(upFlippableIndices.size(), downFlippableIndices.size()); i++) {
            Index upIndex = upFlippableIndices.get(i);
            Index downIndex = downFlippableIndices.get(i);
            if (! upIndex.isNeighbor(downIndex)) {
                double upProp = lattice.flipFaceWeights[upIndex.x][upIndex.y];
                double prop = (faceStates[upIndex.x][upIndex.y] == ns) ? upProp : 1/upProp;
                double downProp = lattice.flipFaceWeights[downIndex.x][downIndex.y];
                prop *= (faceStates[upIndex.x][upIndex.y] == ns) ? downProp : 1/downProp;
                prop *= acceptanceRatioConstant;
                if (rand.nextDouble() < prop) {
                    if (flippableDirection(upIndex) == 1 && flippableDirection(downIndex) == -1) {
                        flipFace(upIndex);
                        flipFace(downIndex);
                    }
                }
            }
        }
    }

    public boolean isDimer(Index coords, int direction) {
        // returns true if there is a dimer in the given direction from the coords face.
        // direction in {0,1,2,3}
        return ((faceStates[coords.x][coords.y] >> direction) & 1) != 0;
    }

    public void flipFace(Index ind) {
        if (faceStates[ind.x][ind.y] == ns || faceStates[ind.x][ind.y] == we) {
            int direction = flippableDirection(ind.x, ind.y);
            currentVolume += 4 * direction;

            faceStates[ind.x][ind.y] ^= 0b1111;
            faceStates[ind.x-1][ind.y] ^= 0b0100;
            faceStates[ind.x][ind.y-1] ^= 0b1000;
            faceStates[ind.x+1][ind.y] ^= 0b0001;
            faceStates[ind.x][ind.y+1] ^= 0b0010;
        }
    }

    public void flipFaceExclusive(int i, int j) {
        if (faceStates[i][j] == ns || faceStates[i][j] == we) {
            int direction = flippableDirection(i, j);
            currentVolume += 4 * direction;

            faceStates[i][j] ^= 0b1111;
        }
    }

    public void consolidateFaceStateFromNeighbors(int i, int j) {
        faceStates[i][j] = 0b0000;
        if(i+1 < lattice.N){
            faceStates[i][j] |= (faceStates[i+1][j] & 0b0001) << 2;
        }
        if(i-1 >= 0) {
            faceStates[i][j] |= (faceStates[i-1][j] & 0b0100) >> 2;
        }
        if(j+1 < lattice.M) {
            faceStates[i][j] |= (faceStates[i][j+1] & 0b0010) << 2;
        }
        if(j-1 >= 0) {
            faceStates[i][j] |= (faceStates[i][j-1] & 0b1000) >> 2;
        }
    }

    public void flipFace(int i, int j) {
        flipFace(new Index(i,j));
    }

    private int flippableDirection(Index ind) {
        return flippableDirection(ind.x, ind.y);
    }

    private int flippableDirection(int i, int j) {
        // returns 1 if flippable in + direction, -1 if in - direction and 0 if not flippable
        int flipDir = 0;
        if (faceStates[i][j] == ns) {
            flipDir = 2 * ((i + j + 1) % 2) - 1;
        } else if (faceStates[i][j] == we) {
            flipDir = 2 * ((i + j) % 2) - 1;
        }
        return flipDir;
    }

    private int computeVolume() {
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
        int reportFreq = numSteps/10;
        long time = System.currentTimeMillis();
        for (int i = 0; i < numSteps; i++) {
            if ((i % reportFreq == 0)) {

                System.out.println("Done with " + i + " steps." + " Average time per 1000 markovSteps: " + (System.currentTimeMillis() - time) * 1000 / reportFreq);
                time = System.currentTimeMillis();
            }
            // Choose parity at random at each step. Can probably just alternate too? 
            // Do we need to take weights into account here?
            markovStep(rand.nextInt(2));
        }
        for (MarkovSimZ2Worker worker : markovWorkers) {
            worker.stopRunning();
        }
    }

    public void simulate(int numSteps, double averageNormalizedHeight) {
        // volumeConstraint is per square. So we multiply it by N*M here.
        averageNormalizedHeight *= lattice.N;
        int volumeGoal = (int) (averageNormalizedHeight * lattice.N * lattice.M);
        int maxSteps = 1000;
        // First get to the desired volume
        for (int i = 0; i < maxSteps; i++) {
            markovStepUpVol(rand.nextInt(2), volumeGoal);
            if (currentVolume == volumeGoal) {
                computeHeightFunctionFromFaceStates();
                break;
            }
        }
        System.out.println("volume " + currentVolume + " reached. Starting simulation.");
        // Now start random updates while maintaining volume.
        int reportFreq = numSteps/10;

        for (int i = 0; i < numSteps; i++) {
            if ((i % reportFreq == 0)) {
                System.out.println("Done with " + i + " steps. Current volume: " + currentVolume);
            }
            markovStepSameVol();
        }
        computeHeightFunctionFromFaceStates();
    }

    public Integer getHeight(int i, int j) {
        if (insideBoundary[i][j]) {
            return heightFunction[i][j];
        } else{
            return 0;
        }
    }

    // Can initialize in many different ways. Should have different versions of this function for different boundary conditions.
    public void initializeFlatSquare() {
        // This corresponds to a matching of only horizontal 
        faceStates = new byte[lattice.N][lattice.M];
        for (int i = 0; i < lattice.N; i++) {
            for (int j = 0; j < lattice.M; j++) {
                // Boundary is boundary of the lattice.
                insideBoundary[i][j] = !((i==0) || (i==lattice.N - 1) || (j==0) || (j==lattice.M - 1));
                if (i % 2 == 1) {
                    // Boundary cases
                    if (j == 0) {
                        faceStates[i][j] = 0b1000;
                    } else if (j ==lattice.M - 1) {
                        faceStates[i][j] = 0b0010;
                    } else if (insideBoundary[i][j]){
                        faceStates[i][j] = ns;
                    }
                }
            }
        }
        computeHeightFunctionFromFaceStates();
        currentVolume = computeVolume();
    }

    public void initializeAztecDiamond() {
        // We assume N = M both odd. Then builds an N-2 x N-2 aztec diamond on the inside.
        faceStates = new byte[lattice.N][lattice.M];
        Index diamondCenter = new Index(lattice.N / 2, lattice.N / 2);
        Index bottomIndex = new Index(lattice.N / 2, 1);
        for (int i = 0; i < lattice.N; i++) {
            for (int j = 0; j < lattice.M; j++) {
                // Boundary of the lattice.
                insideBoundary[i][j] = diamondCenter.l1Dist(i, j) <= lattice.N / 2 - 1;
                if (j < diamondCenter.y) {
                    if (bottomIndex.minus(i, j).isEven()) {
                        faceStates[i][j] = 0b0010;
                    } else {
                        faceStates[i][j] = 0b1000;
                    }
                } else if (j > diamondCenter.y) {
                    if (bottomIndex.minus(i, j).isEven()) {
                        faceStates[i][j] = 0b1000;
                    } else {
                        faceStates[i][j] = 0b0010;
                    }
                } else {
                    if (bottomIndex.minus(i, j).isEven()) {
                        faceStates[i][j] = 0b1010;
                    }
                }
            }
        }
        computeHeightFunctionFromFaceStates();
        currentVolume = computeVolume();
    }



    private void computeHeightFunctionFromFaceStates() {
        // computes the height function associated to the current face states.
        // Note that this only produces an actual height function if the face states encode a perfect matching! - (Could check for this - maybe todo)
        for (int i = 0; i < heightFunction.length; i++) {
            if(i==0) {
                heightFunction[0][0] = 0;
            } else {
                heightFunction[i][0] = heightFunction[i-1][0] + (2 * (i % 2) - 1);
            }
            for (int j = 1; j < heightFunction[i].length; j++) {
                int heightChange = 1 - 4 * ((faceStates[i][j] >> 1) & 1);
                int sign = 2 * ((i+j+1) % 2) - 1;
                heightFunction[i][j] = heightFunction[i][j-1] + sign * heightChange;
            }
        }
    }
}

class Index implements Serializable{
    final int x;
    final int y;

    Index(int x, int y) {this.x=x;this.y=y;}

    public boolean isNeighbor (Index other) {
        return (Math.abs(x - other.x) + Math.abs(y - other.y)) == 1;
    }

    public boolean equals (Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof Index)) {
            return false;
        }
        Index c = (Index) o;
        return (c.x == x) & (c.y == y);
    }

    public Index minus(int i, int j) {
        return new Index(x - i, y - j);
    }

    public int l1Dist(int i, int j) {
        return Math.abs(i - x) + Math.abs(j - y);
    }

    public boolean isEven() {
        return ((x + y) % 2) == 0;
    }
}