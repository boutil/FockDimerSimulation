package dimerSim;

import org.jzy3d.analysis.AWTAbstractAnalysis;
import org.jzy3d.analysis.AnalysisLauncher;
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
import com.jogamp.opengl.awt.GLCanvas;

public class VisualizationZ2 {

    public MarkovSimZ2 sim;


    public VisualizationZ2(MarkovSimZ2 sim) {
        this.sim = sim;
    }

    public void visualizeWeights() {
        // Func3D func = new Func3D((x, y) -> sim.lattice.getUnflippedFaceWeight(x.intValue(), y.intValue()));
        Func3D func = new Func3D((x, y) -> sim.lattice.flipFaceWeights[x.intValue()][y.intValue()]);
        Range range = new Range(0, sim.lattice.N - 1);
        int steps = sim.lattice.N;
    
        // Create the object to represent the function over the given range.
        final Shape surface = new SurfaceBuilder().orthonormal(new OrthonormalGrid(range, steps), func);
    
        surface.setColorMapper(new ColorMapper(new ColorMapRainbow(), surface, new Color(1, 1, 1, .5f)));
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


    public void visualizeSim() {
        Func3D func = new Func3D((x, y) -> sim.getHeight(x.intValue(), y.intValue()).doubleValue());
        Range range = new Range(0, sim.lattice.N - 1);
        int steps = sim.lattice.N;
    
        // Create the object to represent the function over the given range.
        final Shape surface = new SurfaceBuilder().orthonormal(new OrthonormalGrid(range, steps), func);
    
        surface.setColorMapper(new ColorMapper(new ColorMapRainbow(), surface, new Color(1, 1, 1, .5f)));
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


}
