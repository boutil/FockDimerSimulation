package lattices;

import java.io.Serializable;
import java.util.Arrays;

import de.jtem.mfc.field.Complex;

public class HexLattice extends Lattice{
    // Hex coordinates used:

    // Offset coordinates:
    //  x * E + y * NE
    //    / \
    //   |0,1|
    //  / \ / \
    // |0,0|1,0|
    //  \ / \ /
    //   |1,-1|
    //    \ /

    // Cube coordinates:
    // q + r + s = 0

    //      /    \ /    \
    //     |-1,1,0|0,1,-1|
    //  /    \  /   \  /   \
    // |-1,0,1|0,0,0 |1,0,-1|
    //  \    / \    /  \   /
    //     |0,-1,1|1,-1,0|
    //      \    /  \   /

    // Perfect N hexagon is marked by max(|q|,|r|,|s|) <= N
    // Cube => Offset:  (x,y)   = (q,r)
    // Offset ==> Cube: (q,r,s) = (x,y, -x-y)

    // Offset coordinates used for array manipulations
    // Cube coordinates are more natural



    // The flipFaceWeights here are defined as the crossratio of (E * NW * SW)/(NE * W * SE). 
    // !!CAREFUL!! As opposed to the Z2 lattice case, there is only one type of hexagon here so this is also the faceWeight for all faces.
    //  However it allows to do the flip easily in the simulation without considering many cases.

    // (centerX, centerY) corresponds to (0,0,0) in cube coordinates
    private int centerX, centerY;

    public HexLattice(int N, int M, double faceWeightMultiplier) {
        super(N, M, faceWeightMultiplier);

        centerX = N/2;
        centerY = M/2;

        for(double[] row : flipFaceWeights){
            Arrays.fill(row, 1);
        }
        for (int i = 0; i < flipFaceWeights.length; i++) {
            for (int j = 0; j < flipFaceWeights[i].length; j++) {
                flipFaceWeights[i][j] = ((i+j)%2 == 0) ? faceWeightMultiplier : 1/faceWeightMultiplier;
                // flipFaceWeights[i][j] = volumeConstraint;
            }
        }
    }

    public HexLattice(int N, int M) {
        this(N, M, 1);
    }

    public double getUnflippedFaceWeight(int i, int j) {
        return flipFaceWeights[i][j];
    }

    public double getFlipFaceWeightCubeCoords(int q, int r, int s) {
        // q,r,s cube coordinates. return the corresponding flipFaceWeight
        return flipFaceWeights[q + centerX][r + centerY];
    }

    public int[] getCubeCoords(int x, int y) {
        return new int[]{x - centerX, y - centerY, centerX + centerY - x - y};
    }
}
