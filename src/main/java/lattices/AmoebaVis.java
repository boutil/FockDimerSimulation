package lattices;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import il.ac.idc.jdt.DelaunayTriangulation;

public class AmoebaVis extends JPanel{
    
    SchottkyDimers schottkyDimers;

    private BufferedImage amoebaImage;
    private BufferedImage boundaryImage;
    Color[] innerOvalColors = {Color.RED, Color.GREEN, Color.BLUE, Color.YELLOW};
    Color[] ovalColors;

    boolean boundaryPointsComputed = false;
    Complex[][] boundaryPoints;

    boolean colorMapComputed = false;

    boolean isZ2Grid;

    private Complex xCoord = new Complex(1, 0);
    private Complex yCoord = new Complex(Math.cos(Math.PI/3), Math.sin(Math.PI/3));

    private int imageWidth = 1000, imageHeight = 1000;

    public AmoebaVis(SchottkyDimers dimers) {
        Random r = new Random();
        schottkyDimers = dimers;
        amoebaImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_3BYTE_BGR);
        boundaryImage = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_3BYTE_BGR);
        ovalColors = new Color[dimers.numAngles + dimers.getNumGenerators() * 2];
        Arrays.fill(ovalColors, Color.BLACK);
        for (int i = 0; i < ovalColors.length; i++) {
            ovalColors[i] = new Color(r.nextInt(255), r.nextInt(255), r.nextInt(255));
        }
        for (int i = 0; i < dimers.getNumGenerators(); i++) {
            ovalColors[dimers.numAngles + i] = innerOvalColors[i];
        }
        isZ2Grid = dimers.getClass().isAssignableFrom(SchottkyDimersQuad.class);
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
        Graphics2D gAztec = boundaryImage.createGraphics();
        gAmoeba.setColor(Color.WHITE);;
        gAmoeba.fillRect( 0, 0, imageWidth, imageHeight);
        gAztec.setColor(Color.WHITE);;
        gAztec.fillRect( 0, 0, imageWidth, imageHeight);

        ComplexFn amoebaMap = x -> schottkyDimers.amoebaMap(x);
        Complex[][] amoebaPoints = extractOvalPoints(amoebaMap);

        // For now disable until I deal with both Hexagonal and Aztec case. Probably separate classes.

        ComplexFn boundaryMap = x -> schottkyDimers.boundaryCurve(x);
        Complex[][] boundaryPoints = extractOvalPoints(boundaryMap);

        drawPoints(amoebaPoints, gAmoeba);
        drawPoints(boundaryPoints, gAztec);


        gAmoeba.dispose();
        gAztec.dispose();
        repaint();
    }

    interface ComplexFn{
        Complex apply(Complex point) throws Exception;
    }

    private Complex getAztecInTiltedCoords(Complex aztecMap) throws Exception{
        // shrink and rotate.
        // if (Math.abs(aztecMap.re) > 1 || Math.abs(aztecMap.im) > 1 ) {
        //     System.out.println(aztecMap);
        // }
        aztecMap.assignTimes(new Complex(-Math.cos(Math.PI/4), Math.sin(-Math.PI/4)));
        // aztecMap.assignTimes(new Complex(Math.cos(-Math.PI/4), Math.sin(-Math.PI/4)));
        aztecMap.assignDivide(Math.sqrt(2) * 2);
        aztecMap.assignPlus(new Complex(0.5, 0.5));
        aztecMap.assignTimes(imageWidth);
        return aztecMap;
    }

    private Complex getHexagonInTiltedCoordinates(Complex hexMap) throws Exception {
        // Apply hexagonal coordinates point in (x,y) coordinates. Then shring and rotate it to fit simulation pic.
        // if (Math.abs(hexMap.re) > 1 || Math.abs(hexMap.im) > 1 ) {
        //     System.out.println("boundary map larger than 1: " + hexMap);
        // }
        hexMap = xCoord.times(hexMap.re).plus(yCoord.times(hexMap.im));
        hexMap.assignTimes(new Complex(Math.cos(Math.PI/3), Math.sin(Math.PI/3)));
        hexMap.assignTimes((imageHeight - 80) / Math.sqrt(3));
        hexMap.assignPlus(new Complex(imageWidth/2, imageHeight/2));
        return hexMap;
    }

    public void drawBoundaryCurves(Graphics2D g, int imageWidth, int imageHeight) {
        this.imageWidth = imageWidth;
        this.imageHeight = imageHeight;
        if(!boundaryPointsComputed){
            ComplexFn boundaryMap;
            if (isZ2Grid) {
                boundaryMap = x -> getAztecInTiltedCoords(schottkyDimers.boundaryCurve(x));
            } else {
                boundaryMap = x -> getHexagonInTiltedCoordinates(schottkyDimers.boundaryCurve(x));
            }
            boundaryPoints = extractOvalPoints(boundaryMap);
            boundaryPointsComputed = true;
        }
        Color[] whiteColors = new Color[schottkyDimers.numAngles + 1 + schottkyDimers.getNumGenerators()];
        Arrays.fill(whiteColors, Color.WHITE);
        g.setStroke(new BasicStroke(5));
        drawPoints(boundaryPoints, g, whiteColors, false);
    }

    public void drawBoundaryMap(Graphics2D g, int imageWidth, int imageHeight) {
        this.imageWidth = imageWidth;
        this.imageHeight = imageHeight;
        if(!boundaryPointsComputed){
            ComplexFn boundaryMap;
            if (isZ2Grid) {
                boundaryMap = x -> getAztecInTiltedCoords(schottkyDimers.boundaryMap(x));
            } else {
                boundaryMap = x -> getHexagonInTiltedCoordinates(schottkyDimers.boundaryMap(x));
            }
            boundaryPointsComputed = true;
        }
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
        ImageIO.write(boundaryImage, "PNG", new File(filePath));
    }


    private Complex[][] extractOvalPoints(ComplexFn f) {
        int numPointsPerSegment = 300;
        Complex[][] points = schottkyDimers.parametrizeRealOvals(numPointsPerSegment);
        int numSegments = points.length;
        Complex[][] pointsAmoebaMapped = new Complex[numSegments][];
        for (int i = 0; i < numSegments; i++){
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
