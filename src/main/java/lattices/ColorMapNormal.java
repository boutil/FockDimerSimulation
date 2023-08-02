package lattices;

import org.jzy3d.colors.Color;
import org.jzy3d.colors.ColorMapper;
import org.jzy3d.maths.Coord3d;

public class ColorMapNormal extends ColorMapper{

    Coord3d[][] normals;
    double[][] heightFunction;

    public ColorMapNormal(double[][] heightFunction) {
        super();
        this.heightFunction = heightFunction;
        calculateNormals();
    }
    
    // public Color getColor(double x, double y, double z, double zMin, double zMax) {
    //     return colorSphere(normals[(int) x][(int) y]);
    // }

    @Override
    public Color getColor(Coord3d coord) {
        return colorSphere(normals[(int) coord.x][(int) coord.y]);
    }

    private Color colorSphere(Coord3d normal) {
        Coord3d col = normal.getNormalizedTo(0.3f).add(0.5f);
        return new Color(col.x, col.y, col.z);
    }

    private void calculateNormals() {
        normals = new Coord3d[heightFunction.length][heightFunction[0].length];
        int stepSize = 4;
        for (int i = 0; i < heightFunction.length; i++) {
            for (int j = 0; j < heightFunction.length; j++) {
                if ((i <= stepSize) || (j <= stepSize) || (i >= heightFunction.length - stepSize - 1) || (j >= heightFunction[i].length - stepSize - 1)) {
                    normals[i][j] = new Coord3d(0, 0, 1);
                } else {
                    double dx = (heightFunction[i + stepSize][j] - heightFunction[i - stepSize][j]) / (stepSize * 2);
                    double dy = (heightFunction[i][j + stepSize] - heightFunction[i][j - stepSize]) / (stepSize * 2);
                    normals[i][j] = new Coord3d(dx, dy, -1);
                }
                normals[i][j].normalizeTo(1f);
            }
        }
        int a = 3;
    }
}
