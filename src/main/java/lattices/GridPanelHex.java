package lattices;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import de.jtem.mfc.field.Complex;
import dimerSim.Index;
import dimerSim.MarkovSimZ2;

public class GridPanelHex extends JPanel{
    private BufferedImage paintImage;
    private BufferedImage weightsImage;

    private int imageWidth, imageHeight;

    private boolean drawAztecCurves = false;
    private AmoebaVis amoebaVis;

    private MarkovSimZ2 sim;

    private Color[] dimerColors = {Color.BLUE, Color.RED, Color.GREEN};
    // Hexagon diameter is scaling pixels
    public int scaling = 10;

    private Complex xDir = new Complex(2,0);
    private Complex yDir = new Complex(Math.cos(Math.PI / 3), Math.sin(Math.PI / 3));


    public GridPanelHex(MarkovSimZ2 sim) {
        super();
        this.sim = sim;
        imageWidth = sim.lattice.N * scaling;
        imageHeight = sim.lattice.M * scaling;
        paintImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_3BYTE_BGR);
        weightsImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_3BYTE_BGR);
    }

    public GridPanelHex(MarkovSimZ2 sim, AmoebaVis amoebaVis) {
        this(sim);
        drawAztecCurves = true;
        this.amoebaVis = amoebaVis;
    }

    @Override
    public void paintComponent(Graphics g){
        // clear the previous painting
        super.paintComponent(g);
        // cast Graphics to Graphics2D
        g.drawImage(paintImage, 0, 0, null);
    }

    public void updatePaint() {
        Graphics2D g2 = paintImage.createGraphics();
        Graphics2D gWeights = weightsImage.createGraphics();

        double maxWeight = 0;
        for (int i = 0; i < sim.lattice.N; i++) {
            for (int j = 0; j < sim.lattice.M; j++) {
                maxWeight = Math.max(maxWeight, sim.lattice.flipFaceWeights[i][j]);
            }
        }
        gWeights.setColor(Color.WHITE);
        gWeights.drawString("maxWeight: " + maxWeight, 30, 30);

        for (int i = 0; i < sim.lattice.N; i++) {
            for (int j = 0; j < sim.lattice.M; j++) {
                if (sim.insideBoundary[i][j]) {
                    gWeights.setColor(new Color((float) (sim.lattice.flipFaceWeights[i][j]/maxWeight), 0f, 0f));
                    // sim.lattice.getUnflippedFaceWeight(i, j)
                    gWeights.fillRect(i * scaling, j * scaling, 1 * scaling, 1 * scaling);
                    for (int dir = 0; dir < dimerColors.length; dir++) {
                        if (sim.isDimer(new Index(i, j), dir)) {
                            g2.setColor(dimerColors[dir]);
                            g2.fill(getKitePol(i, j, dir));
                            g2.setColor(Color.BLACK);
                            g2.draw(getKitePol(i, j, dir));
                        }
                    }
                }
            }
        }
        g2.dispose();
        repaint();
    }
    
    public void save(String filePath) throws IOException{
        ImageIO.write(paintImage, "PNG", new File(filePath));
    }

    public void saveWeightsPic(String filePath) throws IOException {
        ImageIO.write(weightsImage, "PNG", new File(filePath));
    }

    private Complex get2DDir(int dir) {
        double angle = dir * Math.PI / 3;
        return new Complex(Math.cos(angle), Math.sin(angle));
    }


    private Polygon getKitePol(int i, int j, int dir) {
        Complex centerHex = xDir.times(i).plus(yDir.times(j)).times(scaling);
        Complex centerNextHex = centerHex.plus(get2DDir(dir - 1).times(scaling));
        Complex centerOppositeHex = centerHex.plus(get2DDir(dir).times(scaling));
        Complex centerPrevHex = centerHex.plus(get2DDir(dir + 1).times(scaling));
        return new Polygon(new int[]{(int)centerHex.re, (int)centerNextHex.re, (int)centerOppositeHex.re, (int)centerPrevHex.re}, new int[]{(int)centerHex.im, (int)centerNextHex.im, (int)centerOppositeHex.im, (int)centerPrevHex.im}, 4);
    }
}
