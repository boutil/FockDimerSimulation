package dimerSim;

import java.util.Random;
import java.util.stream.IntStream;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

import lattices.HexLattice;

public class MarkovSimHex extends MarkovSim{
    // states of each face. Encodes a dimer configuration. Can compute a height function from this.
    // 2^0 * E + 2^1 * NE + 2^2 * NW + 2^3 * W + 2^4 * SW + 2^5 * SE
    private byte flippableW = 0b101010;
    private byte flippableE = 0b010101;

    public MarkovSimHex(HexLattice lattice) {
        super(lattice);

        init();

        initializeRegularHexagon();
    }

    public MarkovSimHex(HexLattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        super(lattice, faceStates, insideBoundary);
        init();
    }

    @Override
    protected void init() {
        heightFunction = new int[lattice.N][lattice.M];
        long seed = 42; // for reproducability
        rand = new Random(seed);
        maxParity = 3;
        // for clunky parallelization purposes:
        int numThreads = 20;
        int chunkSize = lattice.N / numThreads;
        markovWorkers = new MarkovSimHexWorker[lattice.N / chunkSize + 1];
        for (int i = 0; i < markovWorkers.length; i++) {
            markovWorkers[i] = new MarkovSimHexWorker(this, IntStream.range(i * chunkSize, Math.min(lattice.N, (i+1) * chunkSize)).toArray());
            markovWorkers[i].start();
        }
    }

    @Override
    public void flipFace(Index ind) {
        if (faceStates[ind.x][ind.y] == flippableE || faceStates[ind.x][ind.y] == flippableW) {
            int direction = flippableDirection(ind.x, ind.y);
            currentVolume += 6 * direction;

            faceStates[ind.x][ind.y] ^= 0b111111;
            faceStates[ind.x-1][ind.y] ^= 0b000001;
            faceStates[ind.x][ind.y-1] ^= 0b000010;
            faceStates[ind.x+1][ind.y-1] ^= 0b000100;
            faceStates[ind.x+1][ind.y] ^= 0b001000;
            faceStates[ind.x][ind.y+1] ^= 0b010000;
            faceStates[ind.x-1][ind.y+1] ^= 0b100000;
        }
    }

    @Override
    public void flipFaceExclusive(int i, int j) {
        if (faceStates[i][j] == flippableE || faceStates[i][j] == flippableW) {
            int direction = flippableDirection(i, j);
            currentVolume += 6 * direction;

            faceStates[i][j] ^= 0b111111;
        }
    }

    @Override
    public void consolidateFaceStateFromNeighbors(int i, int j) {
        faceStates[i][j] = 0b000000;
        if(i+1 < lattice.N){
            faceStates[i][j] |= (faceStates[i+1][j] & 0b001000) >> 3;
        }
        if(i-1 >= 0) {
            faceStates[i][j] |= (faceStates[i-1][j] & 0b000001) << 3;
        }
        if(j+1 < lattice.M) {
            faceStates[i][j] |= (faceStates[i][j+1] & 0b010000) >> 3;
        }
        if(j-1 >= 0) {
            faceStates[i][j] |= (faceStates[i][j-1] & 0b00010) << 3;
        }
        if(i - 1 > 0 && j + 1 < lattice.M) {
            faceStates[i][j] |= (faceStates[i-1][j+1] & 0b10000) >> 3;
        }
        if(i + 1 < lattice.N && j - 1 > 0) {
            faceStates[i][j] |= (faceStates[i+1][j-1] & 0b000100) << 3;
        }
    }

    @Override
    public int flippableDirection(int i, int j) {
        // returns 1 if flippable in + direction, -1 if in - direction and 0 if not flippable
        if (faceStates[i][j] == flippableE) {
            return 1;
        } else if (faceStates[i][j] == flippableW) {
            return -1;
        }
        return 0;
    }

    private void initializeRegularHexagon() {
        
    }
}
