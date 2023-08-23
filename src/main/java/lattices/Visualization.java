package lattices;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;

import org.jzy3d.chart.Chart;
import org.jzy3d.chart.factories.AWTChartFactory;
import org.jzy3d.chart.factories.IChartFactory;
import org.jzy3d.colors.Color;
import org.jzy3d.colors.ColorMapper;
import org.jzy3d.colors.colormaps.ColorMapRainbow;
import org.jzy3d.maths.Range;
import org.jzy3d.plot3d.builder.Func3D;
import org.jzy3d.plot3d.builder.SurfaceBuilder;
import org.jzy3d.plot3d.builder.concrete.OrthonormalGrid;
import org.jzy3d.plot3d.primitives.Shape;
import org.jzy3d.plot3d.rendering.canvas.Quality;


import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.util.Arrays;
import de.jtem.riemann.schottky.SchottkyDimers;
import dimerSim.MarkovSim;
import dimerSim.MarkovSimZ2;

public class Visualization {

    public MarkovSim sim;
    private boolean gridIsZ2;
    private AmoebaVis amoebaVis;
    private GridPanel gridPanel;


    public Visualization(MarkovSim sim, SchottkyDimers schottkyDimers) {
        this.sim = sim;
        if (sim.getClass().isAssignableFrom(MarkovSimZ2.class)) {
            gridIsZ2 = true;
        }
        setVis(new AmoebaVis(schottkyDimers));
    }

    public Visualization(SchottkyDimers schottkyDimers) {
        this.amoebaVis = new AmoebaVis(schottkyDimers);
    }

    public Visualization(MarkovSim sim) {
        this(sim, null);
    }

    public void setSim(MarkovSim sim) {
        this.sim = sim;
        if (sim.getClass().isAssignableFrom(MarkovSimZ2.class)) {
            gridIsZ2 = true;
        }
        setVis(amoebaVis);
    }

    private void setVis(AmoebaVis vis) {
        amoebaVis = vis;
        if (gridIsZ2) {
            gridPanel = new GridPanelZ2(sim, amoebaVis);
        } else {
            gridPanel = new GridPanelHex(sim, amoebaVis);
        }
    }



    public void visualizeWeights() {
        // Func3D func = new Func3D((x, y) -> sim.lattice.getUnflippedFaceWeight(getIndices(x, y, 2)[0], getIndices(x, y, 2)[1]));
        // Func3D func = new Func3D((x, y) -> sim.lattice.flipFaceWeights[x.intValue()][y.intValue()]);
        
        // For only white faces:
        Func3D func = new Func3D((x, y) -> sim.lattice.flipFaceWeights[getIndices(x, y, 1)[0]][getIndices(x, y, 1)[1]]);
        // For only black faces:
        // Func3D func = new Func3D((x, y) -> sim.lattice.flipFaceWeights[x.intValue() - (x.intValue() % 2) + 1][y.intValue() - (y.intValue() % 2) + 1]);
        Range range = new Range(0, sim.lattice.N - 2);
        int steps = sim.lattice.N - 1;
    
        // Create the object to represent the function over the given range.
        final Shape surface = new SurfaceBuilder().orthonormal(new OrthonormalGrid(range, steps), func);
    
        visualizeSurface(surface);
    }

    public void visualizeDiscreteAbelMap() {
        int maxRes = 200;
        Z2LatticeFock fockLattice = (Z2LatticeFock) sim.lattice;
        // Func3D func = new Func3D((x, y) -> fockLattice.discreteAbelMap[x.intValue() - (x.intValue() % 2)][y.intValue() - (y.intValue() % 2)].im[0]);
        Func3D func = new Func3D((x, y) -> fockLattice.discreteAbelMap[x.intValue()][y.intValue()].im[0]);
        // For only black faces:
        // Func3D func = new Func3D((x, y) -> sim.lattice.flipFaceWeights[x.intValue() - (x.intValue() % 2) + 1][y.intValue() - (y.intValue() % 2) + 1]);
        Range range = new Range(0, sim.lattice.N - 1);
        int steps = Math.min(sim.lattice.N, maxRes);
    
        // Create the object to represent the function over the given range.
        final Shape surface = new SurfaceBuilder().orthonormal(new OrthonormalGrid(range, steps), func);
    
        visualizeSurface(surface);
    }

    public void visualizeThetaCrossRatio(double rightStep, double topStep) {
        int res = 200;
        Z2LatticeFock fockLattice = (Z2LatticeFock) sim.lattice;

        ComplexVector topStepV = new ComplexVector(new Complex[] {new Complex(0, topStep)});
        ComplexVector rightStepV = new ComplexVector(new Complex[] {new Complex(0, rightStep)});
        Func3D func = new Func3D((x,y) -> fockLattice.getThetaCrossRatio(new ComplexVector(new Complex[] {new Complex(0, x)}), topStepV, rightStepV).re);
        Range range = new Range(0, 20);
        // Func3D func = new Func3D((x,y) -> fockLattice.getThetaCrossRatio(new ComplexVector(new Complex[] {new Complex(0, x * rightStep / 2)}), topStepV, rightStepV).re);
        // Range range = new Range(0, 100);
        int steps = res;
        final Shape surface = new SurfaceBuilder().orthonormal(new OrthonormalGrid(range, steps), func);
        visualizeSurface(surface);
    }

    private int[] getIndices(Double x, Double y, int onlyEven) {
        // returns indices based on x and y. 
        // OnlyEven 0: Just round down
        // OnlyEven 1: Only even fields
        // OnlyEven 2: Only odd fields
        int[] indices = {x.intValue(), y.intValue()};
        if (onlyEven == 0) {
            return indices;
        } else if (onlyEven == 1) {
            if ((indices[0] + indices[1]) % 2 != 0) {
                indices[0] += 1;
            }
        } else {
            if ((indices[0] + indices[1]) % 2 != 1) {
                indices[0] += 1;
            }
        }
        return indices;
    }


    public void visualizeSim(int smoothingFilterSize) {
        // We can do some subsampling here for better graphics performance (and maybe some smoothing too)
        int maxRes = 200;

        double[][] smoothedHeightFunction = getSmoothedHeightFunction(smoothingFilterSize);

        Func3D func = new Func3D((x, y) -> smoothedHeightFunction[x.intValue()][y.intValue()]);
        Range range = new Range(0, sim.lattice.N - 1);
        int steps = Math.min(sim.lattice.N, maxRes);
    
        // Create the object to represent the function over the given range.
        final Shape surface = new SurfaceBuilder().orthonormal(new OrthonormalGrid(range, steps), func);
        // visualizeSurface(surface);

        surface.setColorMapper(new ColorMapNormal(smoothedHeightFunction));
        surface.setFaceDisplayed(true);
        surface.setWireframeDisplayed(true);
        surface.setWireframeColor(Color.BLACK);
    
        // Create a chart
        IChartFactory f = new AWTChartFactory();
        Chart chart = f.newChart(Quality.Advanced().setHiDPIEnabled(true));
        chart.getScene().getGraph().add(surface);
        chart.open();
        chart.addMouse();

    }


    private void visualizeSurface(final Shape surface) {
        surface.setColorMapper(new ColorMapper(new ColorMapRainbow(), surface, new Color(1, 1, 1, .5f)));
        // surface.setColorMapper(new ColorMapNormal());
        surface.setFaceDisplayed(true);
        surface.setWireframeDisplayed(true);
        surface.setWireframeColor(Color.BLACK);
    
        // Create a chart
        IChartFactory f = new AWTChartFactory();
        Chart chart = f.newChart(Quality.Advanced().setHiDPIEnabled(true));
        chart.getScene().getGraph().add(surface);
        chart.open();
        chart.addMouse();
    }


    public void visualizeDimerConfiguration() {
        // Draw 4 different colored rectangles for the possible different domino types.
        JFrame f = new JFrame();
        Container c = f.getContentPane();

        c.setLayout(new BorderLayout());
        GridPanel p;
        if (gridIsZ2) {
            p = new GridPanelZ2(sim);
        } else {
            p = new GridPanelHex(sim);
        }
        p.updatePaint();
        c.add(p); 
        f.setSize(sim.lattice.N * p.scaling + 100, sim.lattice.M * p.scaling + 100);    
        // make the JFrame visible
        f.setVisible(true);    
        // sets close behavior; EXIT_ON_CLOSE invokes System.exit(0) on closing the JFrame
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);    
    }

    public void saveDimerConfPic(String filePath) throws IOException {
        saveDimerConfPic(filePath, false);
    }

    public void saveDimerConfPic(String filePath, boolean drawCurve) throws IOException {
        gridPanel.drawBoundaryCurves = drawCurve;
        gridPanel.updatePaint();
        gridPanel.save(filePath);
    }

    public void saveDimerConfPic(SchottkyDimers dimers, String filePath) throws IOException {
        setVis(new AmoebaVis(dimers));
        gridPanel.updatePaint();
        gridPanel.save(filePath);
    }

    public void saveAmoebaPic(SchottkyDimers dimers, String filePath) throws IOException {
        setVis(new AmoebaVis(dimers));
        amoebaVis.updatePaint();
        amoebaVis.saveAmoeba(filePath);
    }

    public void saveAztecPic(SchottkyDimers dimers, String filePath) throws IOException {
        setVis(new AmoebaVis(dimers));
        amoebaVis.updatePaint();
        amoebaVis.saveAztec(filePath);
    }

    public void saveWeightsPic(String filepPath) throws IOException {
        GridPanelZ2 p = new GridPanelZ2(sim);
        p.updatePaint();
        p.saveWeightsPic(filepPath);
    }

    public void visualizeAmoeba(SchottkyDimers dimers) {
        JFrame f = new JFrame();
        Container c = f.getContentPane();

        c.setLayout(new BorderLayout());
        AmoebaVis p = new AmoebaVis(dimers);
        p.updatePaint();
        c.add(p);
        f.setSize(1050, 1050);
        // make the JFrame visible
        f.setVisible(true);    
        // sets close behavior; EXIT_ON_CLOSE invokes System.exit(0) on closing the JFrame
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);    
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
