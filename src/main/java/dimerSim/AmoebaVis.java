package dimerSim;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyDimers;

public class AmoebaVis extends JPanel{
    
    SchottkyDimers schottkyDimers;

    private BufferedImage paintImage;
    Color[] ovalColors = {Color.BLACK, Color.BLACK, Color.BLACK, Color.BLACK, Color.RED, Color.GREEN, Color.BLUE};

    private int imageWidth = 1000, imageHeight = 1000;

    public AmoebaVis(SchottkyDimers dimers) {
        schottkyDimers = dimers;
        paintImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_3BYTE_BGR);
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
        g2.setColor(Color.WHITE);;
        g2.fillRect( 0, 0, imageWidth, imageHeight);

        Complex[][] ovalPoints = extractOvalPoints();

        double minRe = Double.MAX_VALUE, maxRe = Double.MIN_VALUE, minIm = Double.MAX_VALUE, maxIm = Double.MIN_VALUE;
        for (int i = 0; i < ovalPoints.length; i++) {
            for (int j = 0; j < ovalPoints[i].length; j++) {
                if(Double.isNaN(ovalPoints[i][j].re) || Double.isNaN(ovalPoints[i][j].im) || Double.isInfinite(ovalPoints[i][j].re) || Double.isInfinite(ovalPoints[i][j].im)) {
                    continue;
                }
                minRe = Math.min(minRe, ovalPoints[i][j].re);
                maxRe = Math.max(maxRe, ovalPoints[i][j].re);
                minIm = Math.min(minIm, ovalPoints[i][j].im);
                maxIm = Math.max(maxIm, ovalPoints[i][j].im);
            }
        }

        for (int i = 0; i < ovalPoints.length; i++) {
            g2.setColor(ovalColors[i]);
            int[] xCoords = new int[ovalPoints[i].length];
            int[] yCoords = new int[ovalPoints[i].length];
            for (int j = 0; j < ovalPoints[i].length; j++) {
                xCoords[j] = (int) (((ovalPoints[i][j].re - minRe) / (maxRe - minRe)) * imageWidth);
                yCoords[j] = (int) (((ovalPoints[i][j].im - minIm) / (maxIm - minIm)) * imageHeight);
            }
            g2.drawPolyline(xCoords, yCoords, ovalPoints[i].length);
        }

        g2.dispose();
        repaint();
    }

    private Complex[][] extractOvalPoints() {
        int numPointsPerSegment = 500;
        Complex[][] points = schottkyDimers.parametrizeRealOvals(numPointsPerSegment);
        int numSegments = schottkyDimers.angles.length + schottkyDimers.getNumGenerators();
        Complex[][] pointsAmoebaMapped = new Complex[numSegments][];
        for (int i = 0; i < points.length; i++){
            List<Complex> mappedP = new LinkedList<Complex>();
            for (int j = 0; j < points[i].length; j++){
                try {
                    Complex amoebaPoint = schottkyDimers.amoebaMap(points[i][j]);
                    if (amoebaPoint.isNaN() || amoebaPoint.isInfinite()) {
                        continue;
                    }
                    mappedP.add(schottkyDimers.amoebaMap(points[i][j]));
                } catch (Exception e) {
                    // System.out.println("While calculating pointsAmoebaMapped: " + e.getMessage() + "P: " + points[i][j]);
                    // if (j != 0) { 
                        //     mappedP[j] = mappedP[j-1];
                    // }
                }
            }
            Complex[] mappedArray = new Complex[mappedP.size()];
            int j = 0;
            for (Complex d : mappedP) {
                mappedArray[j] = d;
                j++;
            }
            pointsAmoebaMapped[i] = mappedArray;
        }
        return pointsAmoebaMapped;
    }

}
