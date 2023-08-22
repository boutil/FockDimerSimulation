package dimerSim;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.stream.IntStream;

import lattices.Z2Lattice;


public class MarkovSimZ2 extends MarkovSim{
    // states of each face. Encodes a dimer configuration. Can compute a height function from this.
    // 2^0 * W + 2^1 * S + 2^2 * E + 2^3 * N
    private byte ns = 0b1010;
    private byte we = 0b0101;

    public MarkovSimZ2(Z2Lattice lattice, boolean flat) {
        super(lattice);

        init();
        
        if (flat){
            initializeFlatSquare();
        } else {
            initializeAztecDiamond();
        }
    }

    public MarkovSimZ2(Z2Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        super(lattice, faceStates, insideBoundary);
        init();
    }

    @Override
    protected void init() {
        heightFunction = new int[lattice.N][lattice.M];
        long seed = 42; // for reproducability
        rand = new Random(seed);
        maxParity = 2;
        // for clunky parallelization purposes:
        int numThreads = 20;
        int chunkSize = lattice.N / numThreads;
        markovWorkers = new MarkovSimZ2Worker[lattice.N / chunkSize + 1];
        for (int i = 0; i < markovWorkers.length; i++) {
            markovWorkers[i] = new MarkovSimZ2Worker(this, IntStream.range(i * chunkSize, Math.min(lattice.N, (i+1) * chunkSize)).toArray());
        }
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

    @Override
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

    @Override
    public void flipFaceExclusive(int i, int j) {
        if (faceStates[i][j] == ns || faceStates[i][j] == we) {
            int direction = flippableDirection(i, j);
            currentVolume += 4 * direction;

            faceStates[i][j] ^= 0b1111;
        }
    }

    @Override
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

    @Override
    public int flippableDirection(int i, int j) {
        // returns 1 if flippable in + direction, -1 if in - direction and 0 if not flippable
        int flipDir = 0;
        if (faceStates[i][j] == ns) {
            flipDir = 2 * ((i + j + 1) % 2) - 1;
        } else if (faceStates[i][j] == we) {
            flipDir = 2 * ((i + j) % 2) - 1;
        }
        return flipDir;
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


    @Override
    protected void computeHeightFunctionFromFaceStates() {
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