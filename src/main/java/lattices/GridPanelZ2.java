package lattices; 

import dimerSim.Index;
import dimerSim.MarkovSim;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.jzy3d.maths.Coord3d;

import de.jtem.numericalMethods.util.Arrays;


public class GridPanelZ2 extends GridPanel {

    // private Color[] dimerColors = {Color.BLUE, Color.RED, Color.GREEN, Color.YELLOW};
    private Color[] dimerColors = {new Color(6, 57, 112), new Color(135,62,35), new Color(118,181,197), new Color(255,200,87)};
    Coord3d[] heightColors = {new Coord3d(6, 57, 112), new Coord3d(135,62,35), new Coord3d(118,181,197), new Coord3d(255,200,87)};

    public int scaling = 4;

    protected BufferedImage heightImage;
    Coord3d[][] normals;

    public GridPanelZ2(MarkovSim sim) {
        super(sim);
        imageWidth = sim.lattice.N * scaling;
        imageHeight = sim.lattice.M * scaling;
        paintImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_4BYTE_ABGR);
        weightsImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_4BYTE_ABGR);
        heightImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_4BYTE_ABGR);
    }

    public GridPanelZ2(MarkovSim sim, AmoebaVis amoebaVis) {
        this(sim);
        drawBoundaryCurves = true;
        this.amoebaVis = amoebaVis;
    }

    public void saveHeightPic(String filePath) throws IOException {
        ImageIO.write(heightImage, "PNG", new File(filePath));
    }

    @Override
    public void updatePaint() {
        Graphics2D g2 = paintImage.createGraphics();
        Graphics2D gWeights = weightsImage.createGraphics();
        Graphics2D gHeight = heightImage.createGraphics();

        g2.setColor(new Color(0, 0, 0, 0));
        g2.fillRect( 0, 0, imageWidth, imageHeight);
        gHeight.setStroke(new BasicStroke(3));
        gHeight.setColor(new Color(0, 0, 0, 0));
        gHeight.fillRect( 0, 0, imageWidth, imageHeight);
        gWeights.setStroke(new BasicStroke(3));
        gWeights.setColor(new Color(0, 0, 0, 0));
        gWeights.fillRect( 0, 0, imageWidth, imageHeight);

        double maxWeight = 0;
        for (int i = 0; i < sim.lattice.N; i++) {
            for (int j = 0; j < sim.lattice.M; j++) {
                maxWeight = Math.max(maxWeight, sim.lattice.flipFaceWeights[i][j]);
            }
        }
        gWeights.setColor(Color.WHITE);
        gWeights.drawString("maxWeight: " + maxWeight, 30, 30);

        calculateNormals(getSmoothedHeightFunction(1));

        for (int i = 0; i < sim.lattice.N; i++) {
            for (int j = 0; j < sim.lattice.M; j++) {
                if (sim.insideBoundary[i][j]) {
                    gWeights.setColor(new Color((float) (sim.lattice.flipFaceWeights[i][j]/maxWeight), 0f, 0f));
                    // sim.lattice.getUnflippedFaceWeight(i, j)
                    gWeights.fillRect(i * scaling, j * scaling, 1 * scaling, 1 * scaling);
                    
                    gHeight.setColor(colorSphere(normals[i][j]));
                    gHeight.fillRect(i * scaling, j * scaling, 1 * scaling, 1 * scaling);
                    
                    for (int dir = 0; dir < dimerColors.length; dir++) {
                        int dimerType = getDimerType(new Index(i, j), dir);
                        if (dimerType != 4) {
                            g2.setColor(dimerColors[dimerType]);
                            g2.fill(getDominoRect(i, j, dir));
                            g2.setColor(Color.BLACK);
                            g2.draw(getDominoRect(i, j, dir));
                        }
                    }
                }
            }
        }
        if (drawBoundaryCurves) {
            amoebaVis.drawBoundaryCurves(g2, imageWidth, imageHeight);
            amoebaVis.drawBoundaryCurves(gHeight, imageWidth, imageHeight);
        }

        g2.dispose();
        repaint();
    }

    private Rectangle getDominoRect(int i, int j, int dir) {
        boolean isHorizontal = (dir % 2) == 1;
        int width = scaling * (isHorizontal ? 2 : 1);
        int height = scaling * (isHorizontal ? 1 : 2);
        int x = scaling * (i - (dir != 2 ? 1 : 0));
        int y = scaling * (j + (dir == 3 ? 1 : 0));
        return new Rectangle(x, y, width, height);
    }

    private int getDimerType(Index coords, int dir) {
        if(sim.isDimer(coords, dir)) {
            if(coords.isEven()) {
                return dir;
            } else {
                return (dir + 2) % 4;
            }
        }
        return 4;
    }

    private Color colorSphere(Coord3d normal) {
        // Map normal to slope in normalized square.
        // Then coordinates give color combination. 
        // Specific for square construction. In general would want to use Newton polygon map.
        Coord3d n = normal.getNormalizedTo(1).mul(1, 1, 0);
        n.rotate(45, new Coord3d(0, 0, 1));
        n.mulSelf((float)Math.sqrt(2));
        n.addSelf(0.5f, 0.5f, 0);
        float[] factors = new float[4];
        Coord3d[] sqVerts = new Coord3d[]{new Coord3d(0, 0, 0), new Coord3d(1, 0, 0), new Coord3d(0, 1, 0), new Coord3d(1, 1, 0)};
        Coord3d finalCol = new Coord3d();
        // factors[0] = 1 - n.x + 1 - n.y;
        // factors[1] = n.x + 1 - n.y;
        // factors[2] = n.x + n.y;
        // factors[3] = 1 - n.x + n.y;
        float facSum = 0;
        for (int i = 0; i < heightColors.length; i++) {
            // factors[i] = Math.max(Math.abs(n.x - sqVerts[i].x), Math.abs(n.y - sqVerts[i].y));
            factors[i] = Math.max(0, Math.min(1, 1 - (float) n.distance(sqVerts[i])));
            facSum += factors[i];
        }
        for (int j = 0; j < heightColors.length; j++) {
            finalCol.addSelf(heightColors[j].mul(factors[j] / facSum));
            
        }
        finalCol.divSelf(255);
        return new Color(finalCol.x, finalCol.y, finalCol.z);
    }

    private void calculateNormals(double[][] heightFunction) {
        normals = new Coord3d[heightFunction.length][heightFunction[0].length];
        int stepSize = 4;
        for (int i = 0; i < heightFunction.length; i++) {
            for (int j = 0; j < heightFunction.length; j++) {
                if ((i <= stepSize) || (j <= stepSize) || (i >= heightFunction.length - stepSize - 1) || (j >= heightFunction[i].length - stepSize - 1)) {
                    normals[i][j] = new Coord3d(0, 0, 1);
                } else {
                    double dx = (heightFunction[i + stepSize][j] - heightFunction[i - stepSize][j]) / (stepSize * 4);
                    double dy = (heightFunction[i][j + stepSize] - heightFunction[i][j - stepSize]) / (stepSize * 4);
                    normals[i][j] = new Coord3d(dx, dy, -1);
                }
                normals[i][j].normalizeTo(1f);
            }
        }
    }
    public double[][] getSmoothedHeightFunction(int filterSize) {
        double[][] smoothedHeightFunction = new double[sim.heightFunction.length][sim.heightFunction.length];
        for (int i = 0; i < smoothedHeightFunction.length; i++) {
            for (int j = 0; j < smoothedHeightFunction.length; j++) {
                smoothedHeightFunction[i][j] = sim.heightFunction[i][j];
            }
        }
        double[][] filter = new double[filterSize][filterSize];
        int n = filterSize / 2;
        for (int i = 0; i < filter.length; i++) {
            Arrays.fill(filter[i], 1/Math.pow(filter.length, 2));
        }
        for (int i = n; i < smoothedHeightFunction.length - n; i++) {
            for (int j = n; j < smoothedHeightFunction[i].length - n; j++) {
                smoothedHeightFunction[i][j] = singlePixelConvolution(sim.heightFunction, filter, i, j);
            }
        }
        return smoothedHeightFunction;
    }

    private double singlePixelConvolution(int[][] pic, double[][] filter, int i, int j) {
        // Assumes a square uneven size filter
        int n = filter.length / 2;
        double val = 0;
        for (int k = -n; k <= n; k++) {
            for (int l = -n; l <= n; l++) {
                val += pic[i + k][j + l] * filter[k + n][l + n];
            }
        }
        return val;
    }

}