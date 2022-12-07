package dimerSim;

import java.io.Serializable;
import java.util.Random;

public class MarkovSimZ2 implements Serializable{
    public Z2Lattice lattice;
    public int[][] heigthFunction;
    // Indicates for each face whether it should be updated.
    public boolean[][] insideBoundary;
    private Random rand;

    public MarkovSimZ2(Z2Lattice lattice) {
        this.lattice = lattice;
        // initialize height function in a sensible way
        heigthFunction = new int[lattice.N][lattice.M];
        insideBoundary = new boolean[lattice.N][lattice.M];
        
        initializeFlatSquare();
        long seed = 42; // for reproducability
        rand = new Random(seed);
    }

    private void markovStep(int parity) {
        // update the height function probabilistically at every even or odd face.
        // recognize if a face is flippable by comparing neighbors. 
        // If so, flip with correct propability.
        
        // For now iterate over the entire lattice. Should probably only consider the flippable faces.
        for (int i = 0; i < lattice.N; i++) {
            for (int j = 0; j < lattice.M; j++) {
                if (!insideBoundary[i][j]) {
                    continue;
                }
                if ((i + j) % 2 == parity) {
                    int faceState = getFaceState(i, j);
                    if (faceState == 5) {
                        double a = (1/lattice.flipFaceWeights[i][j])/2;
                        if (rand.nextDouble() < a){
                            int sign = heigthFunction[i + 1][j] - heigthFunction[i][j]; 
                            heigthFunction[i][j] += sign * 4;
                        }
                    }
                    if (faceState == 10) {
                        double a = lattice.flipFaceWeights[i][j]/2;
                        if (rand.nextDouble() < a){
                            int sign = heigthFunction[i][j + 1] - heigthFunction[i][j];
                            heigthFunction[i][j] += sign * 4;
                        }
                    }
                }
            }    
        }
    }

    public void simulate(int numSteps) {
        for (int i = 0; i < numSteps; i++) {
            markovStep(i % 2);
        }
    }

    private int getFaceState(int i, int j){
        // Returns a representation of which sides of the face have dimers adjacent to it.
        // 2^0 * N + 2^1 * W + 2^2 * S + 2^3 * E. Here the letters are 0 or one depending on if there is a dimer crossing in that direction.
        // i.e. 5 if we have dimers on top+bottom and 10 -> These are the flippable states.

        // all these heightfunction-heightfunction expressions should be either of magnitude 1 or 3. 3 represents a dimer, 1 absence thereof.
        // Therefore take abs() and divide by 2.

        // Let's check for errors first:
        int e = heigthFunction[i + 1][j] - heigthFunction[i][j];
        int w = heigthFunction[i - 1][j] - heigthFunction[i][j];
        int n = heigthFunction[i][j + 1] - heigthFunction[i][j];
        int s = heigthFunction[i][j - 1] - heigthFunction[i][j];
        if ((Math.abs(e) != 1) && (Math.abs(e) != 3)) {
            System.out.println("not a height function.");
        }
        if ((Math.abs(w) != 1) && (Math.abs(w) != 3)) {
            System.out.println("not a height function.");
        }
        if ((Math.abs(n) != 1) && (Math.abs(n) != 3)) {
            System.out.println("not a height function.");
        }
        if ((Math.abs(s) != 1) && (Math.abs(s) != 3)) {
            System.out.println("not a height function.");
        }



        int E = Math.abs(heigthFunction[i + 1][j] - heigthFunction[i][j]) / 2;
        int W = Math.abs(heigthFunction[i - 1][j] - heigthFunction[i][j]) / 2;
        int N = Math.abs(heigthFunction[i][j + 1] - heigthFunction[i][j]) / 2;
        int S = Math.abs(heigthFunction[i][j - 1] - heigthFunction[i][j]) / 2;
        return N + 2 * W + 4 * S + 8 * E;
    }

    // Can initialize in many different ways. Should have different versions of this function for different boundary conditions.
    public void initializeFlatSquare() {
        // This corresponds to a matching of only horizontal 
        for (int i = 0; i < lattice.N; i++) {
            for (int j = 0; j < lattice.M; j++) {
                if (j % 2 == 0){
                    heigthFunction[i][j] = (i%2==0)? 0: -1; 
                } else {
                    heigthFunction[i][j] = (i%2==0)? 1: 2;
                }
                // Boundary is boundary of the lattice.
                insideBoundary[i][j] = !((i==0) || (i==lattice.N - 1) || (j==0) || (j==lattice.M - 1));
            }
        }
    }

    public void initializeAztecDiamond() {
        int diamondCenter = lattice.N/2;
        for (int i = 0; i < lattice.N; i++) {
            for (int j = 0; j < lattice.M; j++) {
                int lineOffSet = 2 * (j - Math.max(0, 2*j - lattice.N));
                if (j % 2 == 0){
                    heigthFunction[i][j] = lineOffSet + i%2;
                } else {
                    heigthFunction[i][j] = lineOffSet + 1 - i%2;
                }
                // Boundary is boundary of the diamond.
                insideBoundary[i][j] = ((Math.abs(diamondCenter - i) + Math.abs(diamondCenter - j)) < lattice.N/2 - 1);
            }
        }
    }
}
