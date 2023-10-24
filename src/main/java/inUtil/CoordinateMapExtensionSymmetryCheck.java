package inUtil;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import javax.swing.JFrame;
import javax.swing.JPanel;

import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.util.Arrays;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import lattices.GridPanel;
import lattices.GridPanelHex;
import lattices.GridPanelZ2;

public class CoordinateMapExtensionSymmetryCheck {
    
    public SchottkyDimersQuad dimers;
    public Complex[] sampledInterval;
    public Complex[] sampledBoundaryCurve;
    public Complex[] sampledAmoeba;
    public Complex[] slopes;
    public int p0Index;

    public int cuspIndex;
    

    public CoordinateMapExtensionSymmetryCheck(SchottkyDimersQuad dimers) {
        this.dimers = dimers;
    }

    public void preCalc(double[] interval) {
        int numPoints = 10000;
        sampledInterval = new Complex[numPoints];
        sampledBoundaryCurve = new Complex[numPoints];
        sampledAmoeba = new Complex[numPoints];
        slopes = new Complex[numPoints];
        for (int i = 0; i < sampledInterval.length; i++) {
            Complex point = new Complex(interval[0] + (i+1) * (interval[1] - interval[0]) / (numPoints + 1), 0);
            sampledInterval[i] = point;
            sampledBoundaryCurve[i] = dimers.boundaryCurve(point);
            sampledAmoeba[i] = dimers.amoebaMap(point);
            slopes[i] = dimers.getDifferentials(point)[1].divide(dimers.getDifferentials(point)[0]).neg();
        }
        // the zero of dxi1.
        double dXi1Min = Double.MAX_VALUE;
        p0Index = 0;
        double Xi1Max = -Double.MAX_VALUE;
        int maxAmoebaIndex = 0;
        double boundaryCurveMax = -Double.MAX_VALUE;
        cuspIndex = 0;
        for (int i = 0; i < numPoints; i++) {
            double dXi1 = dimers.getDifferentials(sampledInterval[i])[0].abs();
            if (dXi1 < dXi1Min) {
                p0Index = i;
                dXi1Min = dXi1;
            }  if (sampledAmoeba[i].re > Xi1Max) {
                Xi1Max = sampledAmoeba[i].re;
                maxAmoebaIndex = i;
            }
            if (sampledBoundaryCurve[i].im > boundaryCurveMax) {
                cuspIndex = i;
                boundaryCurveMax = sampledBoundaryCurve[i].im;
            }
        }
        System.out.println(maxAmoebaIndex + ", " + cuspIndex);
    }

    public Complex[] findPointsAtDistFromAmoebaPoint(int pointIndex, double distance) {
        // Finds points Q1, Q2 that are the intersection of the amoeba with line l parallel to tangent line at sampledAmoeba[pointIndex], where l is distance away from the tangent line.
        double eps = 1E-7;
        
        // Find points with amoebaMap(P) = xi1 on both sides of P0. Do a binary search
        Line l = new Line(sampledAmoeba[pointIndex].minus(distance), -slopes[pointIndex].re);
        Complex[] closestPoints = new Complex[2];
        Complex P1 = sampledInterval[0];
        Complex P2 = sampledInterval[pointIndex];
        double sign = Math.signum(l.dist(dimers.amoebaMap(P1)) - l.dist(dimers.amoebaMap(P2)));
        double currDist = Double.MAX_VALUE;

        while(currDist > eps) {
            Complex Pnew = P1.plus(P2).divide(2);
            if (l.dist(dimers.amoebaMap(Pnew)) * sign < 0) {
                P2 = Pnew;
            } else {
                P1 = Pnew;
            }
            currDist = Math.abs(l.dist(dimers.amoebaMap(Pnew)));
        }
        closestPoints[0] = P1.plus(P2).divide(2);

        P1 = sampledInterval[pointIndex];
        P2 = sampledInterval[sampledInterval.length - 1];
        sign = Math.signum(l.dist(dimers.amoebaMap(P1)) - l.dist(dimers.amoebaMap(P2)));
        currDist = Double.MAX_VALUE;

        while(currDist > eps) {
            Complex Pnew = P1.plus(P2).divide(2);
            if (l.dist(dimers.amoebaMap(Pnew)) * sign < 0) {
                P2 = Pnew;
            } else {
                P1 = Pnew;
            }
            currDist = Math.abs(l.dist(dimers.amoebaMap(Pnew)));
        }
        closestPoints[1] = P1.plus(P2).divide(2);

        return closestPoints;
    }
    
    

    public Complex[] findPointsWithXi1Coord(double xi1, double[] interval) {
        double eps = 1E-7;
        
        // Find points with amoebaMap(P) = xi1 on both sides of P0. Do a binary search
        Complex[] closestPoints = new Complex[2];
        Complex P1 = sampledInterval[0];
        Complex P2 = sampledInterval[p0Index];
        Complex asd = dimers.amoebaMap(P1).minus(dimers.amoebaMap(P2));
        double sign = Math.signum(dimers.amoebaMap(P1).minus(dimers.amoebaMap(P2)).re);
        double currDist = Double.MAX_VALUE;

        while(currDist > eps) {
            Complex Pnew = P1.plus(P2).divide(2);
            double xi1New = dimers.amoebaMap(Pnew).re;
            if ((xi1New - xi1) * sign < 0) {
                P2 = Pnew;
            } else {
                P1 = Pnew;
            }
            currDist = Math.abs(xi1New - xi1);
        }
        closestPoints[0] = P1.plus(P2).divide(2);

        P1 = sampledInterval[p0Index];
        P2 = sampledInterval[sampledInterval.length - 1];
        sign = Math.signum(dimers.amoebaMap(P1).minus(dimers.amoebaMap(P2)).re);
        currDist = Double.MAX_VALUE;

        while(currDist > eps) {
            Complex Pnew = P1.plus(P2).divide(2);
            double xi1New = dimers.amoebaMap(Pnew).re;
            if ((xi1New - xi1) * sign < 0) {
                P2 = Pnew;
            } else {
                P1 = Pnew;
            }
            currDist = Math.abs(xi1New - xi1);
        }
        closestPoints[1] = P1.plus(P2).divide(2);

        return closestPoints;
    }

    // Given P1, P2 on real oval, find their parallel lines in (x,y) plane and compute their intersection.
    public Complex findParallelsIntersection(Complex P1, Complex P2) {
        Complex Q1 = dimers.boundaryCurve(P1);
        Complex Q2 = dimers.boundaryCurve(P2);

        Complex slope1 = dimers.getDifferentials(P1)[1].divide(dimers.getDifferentials(P1)[0]).neg();
        Complex slope2 = dimers.getDifferentials(P2)[1].divide(dimers.getDifferentials(P2)[0]).neg();

        double s1 = slope1.re;
        double s2 = slope2.re;

        double intersectionX = -Q1.im + Q2.im + (Q1.re * s1) - (Q2.re * s2);
        intersectionX /= s1 - s2;
        double intersectionY = (intersectionX - Q1.re) * s1 + Q1.im;
        
        return new Complex(intersectionX, intersectionY);
    }
    
    
    public static void main(String[] args) {
        // Construct some Harnack data
        // G1LargeHole
        double[] schottkyParamsCol = {0.9, 1, 0.9, -1, 0.00000001, 0};
        double[][] angles = {{-2.4,-2.4}, {-1.4, -1.4}, {-1.3, 2.3}, {2.4, 2.4}};

        // G1LargeHole2Angles
        // double[] schottkyParamsCol = {0.9, 1, 0.9, -1, 0.15, 0};
        // double[][] angles = {{-2.4, -0.5}, {-0.4, -0.4}, {0.5, 1.3}, {1.4, 1.4}};

        SchottkyDimersQuad schottkyDimers = new SchottkyDimersQuad(new SchottkyData(schottkyParamsCol), angles);

        CoordinateMapExtensionSymmetryCheck check = new CoordinateMapExtensionSymmetryCheck(schottkyDimers);

        // Find points with same Xi1 coordinate
        double[] interval = angles[2];
        check.preCalc(interval);

        int numCoords = 3000;
        double[] dists = new double[numCoords];
        double intervalLength = 4;
        for (int i = 0; i < dists.length; i++) {
            // double min = check.dimers.amoebaMap(check.sampledInterval[0]).re + 0.5, max = check.dimers.amoebaMap(check.sampledInterval[check.p0Index]).re;
            dists[i] = i * intervalLength/numCoords;
        }

        Complex[] intersectionPoints = new Complex[dists.length];
        int pointIndex = check.cuspIndex;
        for (int i = 0; i < dists.length; i++) {
            Complex[] points = check.findPointsAtDistFromAmoebaPoint(pointIndex, dists[i]);
    
            // Calculate the intersection of their parallel lines in the (x,y) plane.
            Complex intersectionPoint = check.findParallelsIntersection(points[0], points[1]);
            
            // System.out.println(intersectionPoint);
            intersectionPoints[i] = intersectionPoint;
        }

        visualizeIntersectionCurve(intersectionPoints, check.sampledBoundaryCurve, check.sampledBoundaryCurve[pointIndex]);


        // Check that they all have the same x-coordinate.
    }

    public static void visualizeIntersectionCurve(Complex[] points, Complex[] boundaryPoints, Complex cusp) {
        JFrame f = new JFrame();
        Container c = f.getContentPane();

        c.setLayout(new BorderLayout());

        int frameSize = 1000;
        int[][] x = new int[3][];
        int[][] y = new int[3][];
        x[0] = new int[points.length];
        y[0] = new int[points.length];
        x[1] = new int[boundaryPoints.length];
        y[1] = new int[boundaryPoints.length];
        x[2] = new int[1];
        y[2] = new int[1];

        double minX = Double.MAX_VALUE, maxX = -Double.MAX_VALUE, minY = Double.MAX_VALUE, maxY = -Double.MAX_VALUE;
        for (int i = 0; i < boundaryPoints.length; i++) {
            minX = Math.min(minX, boundaryPoints[i].re);
            maxX = Math.max(maxX, boundaryPoints[i].re);
            minY = Math.min(minY, boundaryPoints[i].im);
            maxY = Math.max(maxY, boundaryPoints[i].im);
        }
        double centerX = (minX + maxX) / 2;
        double centerY = (minY + maxY) / 2;
        double scaling = 2 / Math.max(maxX - minX, maxY - minY);

        for (int i = 0; i < points.length; i++) {
            Complex scaledPoint = points[i].minus(new Complex(centerX, centerY)).times(scaling).conjugate().plus(new Complex(1, 1)).times(frameSize / 2);
            x[0][i] = (int)scaledPoint.re;
            y[0][i] = (int)scaledPoint.im;
        }
        for (int i = 0; i < boundaryPoints.length; i++) {
            Complex scaledBoundaryPoint = boundaryPoints[i].minus(new Complex(centerX, centerY)).times(scaling).conjugate().plus(new Complex(1, 1)).times(frameSize / 2);
            x[1][i] = (int)scaledBoundaryPoint.re;
            y[1][i] = (int)scaledBoundaryPoint.im;
        }
        Complex scaledCusp = cusp.minus(new Complex(centerX, centerY)).times(scaling).conjugate().plus(new Complex(1, 1)).times(frameSize / 2);
        x[2][0] = (int) scaledCusp.re;
        y[2][0] = (int) scaledCusp.im;

        
        LinePanel p = new LinePanel(x, y);

        c.add(p);

        f.setSize(frameSize, frameSize);    
        // make the JFrame visible
        f.setVisible(true);    
        // sets close behavior; EXIT_ON_CLOSE invokes System.exit(0) on closing the JFrame
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
    }

    public static class LinePanel extends JPanel {
        int[][] x;
        int[][] y;

        public LinePanel(int[][] x, int[][] y) {
            super();
            this.x = x;
            this.y = y;
        }

        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D)g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                RenderingHints.VALUE_ANTIALIAS_ON);
            // for (int i = 0; i < x.length; i++) {
            //     g2.drawPolyline(x[i], y[i], x[i].length);
            // }
            g2.drawPolyline(x[1], y[1], x[1].length);
            g2.setColor(Color.RED);
            g2.drawPolyline(x[0], y[0], x[0].length);
            g2.setColor(Color.GREEN);
            g2.drawRect(x[2][0] - 5, y[2][0] - 5, 10, 10);

        }
    }

    public class Line {
        public double a, b, c;

        public Line(Complex point, double slope) {
            a = slope;
            b = -1;
            c = point.im - point.re * slope;
        }

        public double dist(Complex p) {
            // returns signed distance.
            return (a * p.re + b * p.im + c) / Math.sqrt(a * a + b * b);
        }
    }
}


