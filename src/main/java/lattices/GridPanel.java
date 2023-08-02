package lattices;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import org.apache.commons.lang3.NotImplementedException;

import de.jtem.mfc.field.Complex;
import dimerSim.Index;
import dimerSim.MarkovSim;

public class GridPanel extends JPanel{
    protected BufferedImage paintImage;
    protected BufferedImage weightsImage;

    protected int imageWidth, imageHeight;

    protected boolean drawAztecCurves = false;
    protected AmoebaVis amoebaVis;

    protected MarkovSim sim;

    protected Color[] dimerColors = {Color.BLUE, Color.RED, Color.GREEN, Color.YELLOW};
    // Hexagon diameter is scaling pixels
    public int scaling = 4;


    public GridPanel(MarkovSim sim) {
        super();
        this.sim = sim;
    }

    public GridPanel(MarkovSim sim, AmoebaVis amoebaVis) {
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
        throw new NotImplementedException();
    }
    
    public void save(String filePath) throws IOException{
        ImageIO.write(paintImage, "PNG", new File(filePath));
    }

    public void saveWeightsPic(String filePath) throws IOException {
        ImageIO.write(weightsImage, "PNG", new File(filePath));
    }
}
