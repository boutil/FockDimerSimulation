package dimerSim;

public class DimerConfigZ2 {
    
    // This class maintains both a 
    private int[][] heightFunction;
    private boolean[][] insideBoundary;
    private int[][] faceStates;

    // points on boundary ordered counterclockwise.
    private int[][] boundary;

    public DimerConfigZ2(int N, int M){
        heightFunction = new int[N][M];
        insideBoundary = new boolean[N][M];
        faceStates = new int[N][M];
        
    }

    public void flipFace(int i, int j) {
        // Checks if the face is in a flippable state and flips it if possible. Exception otherwise.

    }

}
