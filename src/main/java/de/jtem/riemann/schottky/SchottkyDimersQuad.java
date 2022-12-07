package de.jtem.riemann.schottky;

public class SchottkyDimersQuad extends SchottkyDimers {
    public double alphaMinus;
    public double betaMinus;
    public double alphaPlus;
    public double betaPlus;

    public SchottkyDimersQuad(SchottkyData data, double[] angles) {
        super(data, angles);
        alphaMinus = angles[0]; 
        betaMinus = angles[1];
        alphaPlus = angles[2];
        betaPlus = angles[3];
    }
}
