package inUtil;

import de.jtem.mfc.field.Complex;

public class hexagonalCaseCheck {
    public static void main(String[] args) {
        double angleP = 0.235;
        Complex P = new Complex(Math.cos(angleP), Math.sin(angleP)).times(-1);
        // double[] argAngles = {0, Math.PI / 3 , 2 * Math.PI / 3};
        double[] argAngles = {0, 0.25, 0.5};
        Complex[] angles = {new Complex(Math.cos(argAngles[0]), Math.sin(argAngles[0])), new Complex(Math.cos(argAngles[1]), Math.sin(argAngles[1])), new Complex(Math.cos(argAngles[2]), Math.sin(argAngles[2]))};
        
        int numPoints = 10;
        for(int i = 0; i < numPoints; i++) {
            double currentAngle = i * Math.PI / numPoints;
            Complex currentP = new Complex(Math.cos(currentAngle), Math.sin(currentAngle));

            double[] xyz = getXYZ(currentP, angles);
            double x = xyz[0], y = xyz[1], z = xyz[2];
            // System.out.println(x + ", " + y + ", " + z);
            System.out.println(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2) - x * y - y * z - z * x);
        }

        
    }

    public static double[] getXYZ(Complex P, Complex[] angles) {
        Complex alpha = angles[0];
        Complex beta = angles[1];
        Complex gamma = angles[2];
        Complex dXi1 = P.minus(alpha).invert().plus(P.plus(alpha).invert());
        Complex dXi2 = P.minus(beta).invert().plus(P.plus(beta).invert());
        Complex dXi3 = P.minus(gamma).invert().plus(P.plus(gamma).invert());
        Complex dXi = P.minus(alpha).invert().minus(P.plus(alpha).invert()).minus(P.minus(beta).invert()).plus(P.plus(beta).invert()).plus(P.minus(gamma).invert()).minus(P.plus(gamma).invert());
        // Complex dXi1 = P.minus(alpha).invert().plus(P.plus(alpha).invert());
        // Complex dXi2 = P.minus(beta).invert().plus(P.plus(beta).invert());
        // Complex dXi3 = P.minus(gamma).invert().plus(P.plus(gamma).invert());
        // Complex dXi = P.minus(alpha).invert().minus(P.plus(alpha).invert()).plus(P.minus(beta).invert().minus(P.plus(beta).invert())).plus(P.minus(gamma).invert().minus(P.plus(gamma).invert()));

        Complex R1 = dXi1.minus(dXi3).divide(dXi);
        Complex R2 = dXi2.minus(dXi3).divide(dXi);

        double x = - R2.im / (R1.re * R2.im - R2.re * R1.im);
        double y = R1.im / (R1.re * R2.im - R2.re * R1.im);
        // double x = 1/(- R2.re + R1.re * (R2.im/R1.im));
        // double y = 1/(R1.re - R2.re * (R1.im/R2.im));
        double z = - x - y;

        Complex dF = dXi1.times(x).plus(dXi2.times(y)).plus(dXi3.times(z)).plus(dXi);
        // System.out.println(dF.abs());
        return new double[]{x,y,z};
    }
}
