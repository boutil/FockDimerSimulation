package lattices;

import org.jzy3d.colors.Color;
import org.jzy3d.colors.ColorMapper;
import org.jzy3d.maths.Coord3d;

public class ColorMapNormal extends ColorMapper{

    Coord3d[][] normals;
    double[][] heightFunction;

    // Colors for the 4 principal sides corresponding to frozen regions. Rest will be a mix of them.
    Coord3d[] mainColors = {new Coord3d(6, 57, 112), new Coord3d(135,62,35), new Coord3d(118,181,197), new Coord3d(255,200,87)};

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
        // Map normal to slope in normalized square.
        // Then coordinates give color combination. 
        // Specific for square construction. In general would want to use Newton polygon map.
        Coord3d n = normal.getNormalizedTo(1).mul(1, 1, 0);
        n.rotate(45, new Coord3d(0, 0, 1));
        n.addSelf(0.5f, 0.5f, 0);
        float[] factors = new float[4];
        Coord3d[] sqVerts = new Coord3d[]{new Coord3d(0, 0, 0), new Coord3d(1, 0, 0), new Coord3d(0, 1, 0), new Coord3d(1, 1, 0)};
        Coord3d finalCol = new Coord3d();
        // factors[0] = 1 - n.x + 1 - n.y;
        // factors[1] = n.x + 1 - n.y;
        // factors[2] = n.x + n.y;
        // factors[3] = 1 - n.x + n.y;
        float facSum = 0;
        for (int i = 0; i < mainColors.length; i++) {
            // factors[i] = Math.max(Math.abs(n.x - sqVerts[i].x), Math.abs(n.y - sqVerts[i].y));
            factors[i] = Math.max(0, Math.min(1, 1 - (float) n.distance(sqVerts[i])));
            facSum += factors[i];
        }
        for (int j = 0; j < mainColors.length; j++) {
            finalCol.addSelf(mainColors[j].mul(factors[j] / facSum));
            
        }
        finalCol.divSelf(255);
        return new Color(finalCol.x, finalCol.y, finalCol.z);
    }

    private void calculateNormals() {
        normals = new Coord3d[heightFunction.length][heightFunction[0].length];
        int stepSize = 4;
        for (int i = 0; i < heightFunction.length; i++) {
            for (int j = 0; j < heightFunction.length; j++) {
                if ((i <= stepSize) || (j <= stepSize) || (i >= heightFunction.length - stepSize - 1) || (j >= heightFunction[i].length - stepSize - 1)) {
                    normals[i][j] = new Coord3d(0, 0, 1);
                } else {
                    double dx = (heightFunction[i + stepSize][j] - heightFunction[i - stepSize][j]) / (stepSize * 4);
                    double dy = (heightFunction[i][j + stepSize] - heightFunction[i][j - stepSize]) / (stepSize * 4);
                    normals[i][j] = new Coord3d(dx, dy, -1);
                }
                normals[i][j].normalizeTo(1f);
            }
        }
    }
}
