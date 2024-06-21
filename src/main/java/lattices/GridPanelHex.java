package lattices;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.image.BufferedImage;

import de.jtem.mfc.field.Complex;
import dimerSim.Index;
import dimerSim.MarkovSim;

public class GridPanelHex extends GridPanel{
    // private Color[] dimerColors = {Color.BLUE, Color.RED, Color.GREEN};
    private Color[] dimerColors = {new Color(6, 57, 112), new Color(135,62,35), new Color(118,181,197)};
    // Hexagon diameter is scaling pixels

    private Complex xDir = new Complex(1,0);
    private Complex yDir = new Complex(Math.cos(Math.PI / 3), Math.sin(Math.PI / 3));


    public GridPanelHex(MarkovSim sim) {
        super(sim);
        scaling = 12;
        imageWidth = sim.lattice.N * scaling + sim.lattice.M * scaling / 2;
        imageHeight = (int) (sim.lattice.M * scaling * Math.sqrt(3) / 2);
        paintImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_4BYTE_ABGR);
        weightsImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_4BYTE_ABGR);
    }

    public GridPanelHex(MarkovSim sim, AmoebaVis amoebaVis) {
        this(sim);
        drawBoundaryCurves = true;
        this.amoebaVis = amoebaVis;
    }

    @Override
    public void updatePaint() {
        Graphics2D g2 = paintImage.createGraphics();
        Graphics2D gWeights = weightsImage.createGraphics();

        g2.setColor(new Color(0, 0, 0, 0));
        g2.fillRect( 0, 0, imageWidth, imageHeight);

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
                    // gWeights.setColor(new Color((float) (sim.lattice.flipFaceWeights[i][j]/maxWeight), 0f, 0f));
                    // sim.lattice.getUnflippedFaceWeight(i, j)
                    // gWeights.fillRect(i * scaling, j * scaling, 1 * scaling, 1 * scaling);
                    for (int dir = 0; dir < 6; dir++) {
                        if (sim.isDimer(new Index(i, j), dir)) {
                            g2.setColor(dimerColors[dir % 3]);
                            g2.fill(getKitePol(i, j, dir));
                            g2.setColor(Color.BLACK);
                            g2.draw(getKitePol(i, j, dir));
                        }
                    }
                }
            }
        }
        if (drawBoundaryCurves) {
            amoebaVis.drawBoundaryCurves(g2, imageWidth, imageHeight);
        }
        g2.dispose();
        repaint();
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
