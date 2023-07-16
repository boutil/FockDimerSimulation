package dimerSim;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyDimers;

public class AmoebaVis extends JPanel{
    
    SchottkyDimers schottkyDimers;

    private BufferedImage amoebaImage;
    private BufferedImage aztecImage;
    Color[] innerOvalColors = {Color.RED, Color.GREEN, Color.BLUE};
    Color[] ovalColors;

    private int imageWidth = 1000, imageHeight = 1000;

    public AmoebaVis(SchottkyDimers dimers) {
        schottkyDimers = dimers;
        amoebaImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_3BYTE_BGR);
        aztecImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_3BYTE_BGR);
        ovalColors = new Color[dimers.numAngles + dimers.getNumGenerators()];
        Arrays.fill(ovalColors, Color.BLACK);
        for (int i = 0; i < dimers.getNumGenerators(); i++) {
            ovalColors[dimers.numAngles + i] = innerOvalColors[i];
        }
    }

    @Override
    public void paintComponent(Graphics g){
        // clear the previous painting
        super.paintComponent(g);
        // cast Graphics to Graphics2D
        g.drawImage(amoebaImage, 0, 0, null);
    }

    public void updatePaint() {
        Graphics2D gAmoeba = amoebaImage.createGraphics();
        Graphics2D gAztec = aztecImage.createGraphics();
        gAmoeba.setColor(Color.WHITE);;
        gAmoeba.fillRect( 0, 0, imageWidth, imageHeight);
        gAztec.setColor(Color.WHITE);;
        gAztec.fillRect( 0, 0, imageWidth, imageHeight);

        ComplexFn amoebaMap = x -> schottkyDimers.amoebaMap(x);
        Complex[][] amoebaPoints = extractOvalPoints(amoebaMap);
        ComplexFn aztecMap = x -> schottkyDimers.aztecArcticCurve(x);
        Complex[][] aztecPoints = extractOvalPoints(aztecMap);

        drawPoints(amoebaPoints, gAmoeba);
        drawPoints(aztecPoints, gAztec);


        gAmoeba.dispose();
        gAztec.dispose();
        repaint();
    }

    interface ComplexFn{
        Complex apply(Complex point) throws Exception;
    }

    private Complex getAztecInTiltedCoords(Complex point) throws Exception{
        // Apply aztec map from paper and shrink and rotate it.
        Complex aztecMap = schottkyDimers.aztecArcticCurve(point);
        if (Math.abs(aztecMap.re) > 1 || Math.abs(aztecMap.im) > 1 ) {
            System.out.println(aztecMap);
        }
        aztecMap.assignTimes(new Complex(Math.cos(Math.PI/4), Math.sin(Math.PI/4)));
        // aztecMap.assignTimes(new Complex(Math.cos(-Math.PI/4), Math.sin(-Math.PI/4)));
        aztecMap.assignDivide(Math.sqrt(2) * 2);
        aztecMap.assignPlus(new Complex(0.5, 0.5));
        aztecMap.assignTimes(imageWidth);
        return aztecMap;
    }

    public void drawAztecCurves(Graphics2D g, int imageWidth, int imageHeight) {
        this.imageWidth = imageWidth;
        this.imageHeight = imageHeight;
        ComplexFn aztecMap = x -> getAztecInTiltedCoords(x);
        Complex[][] aztecPoints = extractOvalPoints(aztecMap);
        Color[] whiteColors = new Color[schottkyDimers.numAngles + 1 + schottkyDimers.getNumGenerators()];
        Arrays.fill(whiteColors, Color.WHITE);
        g.setStroke(new BasicStroke(6));
        drawPoints(aztecPoints, g, whiteColors, false);
    }

    private void drawPoints(Complex[][] points, Graphics2D g) {
        drawPoints(points, g, ovalColors, true);
    }

    private void drawPoints(Complex[][] points, Graphics2D g, Color[] ovalColors, boolean normalize) {
        double minRe = Double.MAX_VALUE, maxRe = Double.MIN_VALUE, minIm = Double.MAX_VALUE, maxIm = Double.MIN_VALUE;
        for (int i = 0; i < points.length; i++) {
            for (int j = 0; j < points[i].length; j++) {
                if(Double.isNaN(points[i][j].re) || Double.isNaN(points[i][j].im) || Double.isInfinite(points[i][j].re) || Double.isInfinite(points[i][j].im)) {
                    continue;
                }
                minRe = Math.min(minRe, points[i][j].re);
                maxRe = Math.max(maxRe, points[i][j].re);
                minIm = Math.min(minIm, points[i][j].im);
                maxIm = Math.max(maxIm, points[i][j].im);
            }
        }

        for (int i = 0; i < points.length; i++) {
            g.setColor(ovalColors[i]);
            int[] xCoords = new int[points[i].length];
            int[] yCoords = new int[points[i].length];
            for (int j = 0; j < points[i].length; j++) {
                if (normalize) {
                    xCoords[j] = (int) (((points[i][j].re - minRe) / (maxRe - minRe)) * imageWidth);
                    yCoords[j] = (int) (((points[i][j].im - minIm) / (maxIm - minIm)) * imageHeight);
                } else {
                    xCoords[j] = (int) points[i][j].re;
                    yCoords[j] = (int) points[i][j].im;
                }
            }
            g.drawPolyline(xCoords, yCoords, points[i].length);
        }
    }

    public void saveAmoeba(String filePath) throws IOException{
        ImageIO.write(amoebaImage, "PNG", new File(filePath));
    }

    public void saveAztec(String filePath) throws IOException{
        ImageIO.write(aztecImage, "PNG", new File(filePath));
    }


    private Complex[][] extractOvalPoints(ComplexFn f) {
        int numPointsPerSegment = 500;
        Complex[][] points = schottkyDimers.parametrizeRealOvals(numPointsPerSegment);
        int numSegments = schottkyDimers.numAngles + schottkyDimers.getNumGenerators();
        Complex[][] pointsAmoebaMapped = new Complex[numSegments][];
        for (int i = 0; i < points.length; i++){
            List<Complex> mappedP = new LinkedList<Complex>();
            for (int j = 0; j < points[i].length; j++){
                try {
                    Complex amoebaPoint = f.apply(points[i][j]);
                    if (amoebaPoint.isNaN() || amoebaPoint.isInfinite()) {
                        continue;
                    }
                    mappedP.add(amoebaPoint);
                } catch (Exception e) {
                    System.out.println("While calculating pointsAmoebaMapped: " + e.getMessage() + "P: " + points[i][j]);

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
