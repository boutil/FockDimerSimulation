package dimerSim;

import java.util.Random;
import java.util.stream.IntStream;

import lattices.HexLattice;

public class MarkovSimHex extends MarkovSim{
    // states of each face. Encodes a dimer configuration. Can compute a height function from this.
    // 2^0 * E + 2^1 * NE + 2^2 * NW + 2^3 * W + 2^4 * SW + 2^5 * SE
    private byte flippableW = 0b101010;
    private byte flippableE = 0b010101;

    public double sizeQ, sizeR, sizeS = 1;


    public MarkovSimHex(HexLattice lattice) {
        this(lattice, new double[]{1,1,1});
    }
    
    public MarkovSimHex(HexLattice lattice, double[] sizes) {
        super(lattice);
        this.sizeQ = sizes[0];
        this.sizeR = sizes[1];
        this.sizeS = sizes[2];
        int max = lattice.N / 2 - 4;
        initializeHexagon((int)(max * sizeQ), (int)(max * sizeR), (int)(max * sizeS));
        init();
    }

    public MarkovSimHex(HexLattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        super(lattice, faceStates, insideBoundary);
        init();
    }

    @Override
    protected void init() {
        heightFunction = new int[lattice.N][lattice.M];
        long seed = 42; // for reproducability
        rand = new Random();
        maxParity = 3;
        gpuSim = new GPUSimHex(lattice, faceStates, insideBoundary);
    }
    
    @Override
    protected void createWorkers(int numThreads) {
        int chunkSize = lattice.N / numThreads;
        markovWorkers = new MarkovSimHexWorker[lattice.N / chunkSize + 1];
        for (int i = 0; i < markovWorkers.length; i++) {
            markovWorkers[i] = new MarkovSimHexWorker(this, IntStream.range(i * chunkSize, Math.min(lattice.N, (i+1) * chunkSize)).toArray());
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
            faceStates[i][j] |= (faceStates[i][j-1] & 0b000010) << 3;
        }
        if(i - 1 > 0 && j + 1 < lattice.M) {
            faceStates[i][j] |= (faceStates[i-1][j+1] & 0b100000) >> 3;
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

    private void initializeHexagon(int sizeQ, int sizeR, int sizeS) {
        // Point where the initialization lines meet.
        int[] centralCubeCoords = {sizeR - sizeS, sizeQ - sizeR, sizeS - sizeQ};
        HexLattice lat = (HexLattice) lattice;
        for (int i = 0; i < lattice.N; i++) {
            for (int j = 0; j < lattice.M; j++) {
                int[] qrs = lat.getCubeCoords(i, j);
                boolean inside = (Math.abs(qrs[0]) < sizeQ && Math.abs(qrs[1]) < sizeS && Math.abs(qrs[2]) < sizeR);
                insideBoundary[i][j] = inside;
                // Lots of cases to consider. This is a separation into three regions and the resulting faceStates
                // Right quadrant first
                // Boundaries:
                int[] coordDiff = {qrs[0] - centralCubeCoords[0], qrs[1] - centralCubeCoords[1], qrs[2] - centralCubeCoords[2]};
                // Eastern segment
                if ((coordDiff[0] == 0 && coordDiff[1] >= 0) || (coordDiff[2] == 0 && coordDiff[0] >= 0)) {
                    faceStates[i][j] |= 0b000001;
                }
                // SW segment
                if ((coordDiff[2] == 0 && coordDiff[0] >= 0) || (coordDiff[1] == 0 && coordDiff[2] >= 0)) {
                    faceStates[i][j] |= 0b010000;
                }
                // NW segment
                if ((coordDiff[1] == 0 && coordDiff[2] >= 0) || (coordDiff[0] == 0 && coordDiff[1] >= 0)) {
                    faceStates[i][j] |= 0b000100;
                }
                // Non-boundaries:
                // eastern segment
                if (coordDiff[0] > 0 && coordDiff[2] < 0) {
                    faceStates[i][j] |= 0b001001;
                }
                // NW segment
                if (coordDiff[0] < 0 && coordDiff[1] > 0) {
                    faceStates[i][j] |= 0b100100;
                }
                // SW segment
                if (coordDiff[1] < 0 && coordDiff[2] > 0) {
                    faceStates[i][j] |= 0b010010;
                }
            }
        }
    }

    private void initializeRegularHexagon() {
        int size = lattice.N / 2 - 2;
        initializeHexagon(size, size, size);
    }

    public void toPlanePartition() {
        // transforms this configuration into a plane partition
        // We normalize our height function in the following way:
        // planePartition[0][0] corresponds to the bottom left of the hexagon
        // Edges of type NE, SW do not change the height.
        // Edges of type E, W increase height by one going from left to right
        // Edges of type NW, SE increase height by one going in NW direction.
        int[][] planePartition;

    }

}
