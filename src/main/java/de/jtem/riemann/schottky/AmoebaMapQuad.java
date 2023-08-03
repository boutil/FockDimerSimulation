package de.jtem.riemann.schottky;

import de.jtem.mfc.field.Complex;

public class AmoebaMapQuad extends AmoebaMap{

    AmoebaMapQuad(SchottkyDimers schottky, Complex P0) {
        super(schottky, P0);
    }

    @Override
    protected void calculateIncrements(final SchottkyGroupElement element) {
        for (Complex angle : schottky.getAngles()[0]) {
            Complex diff = sP.minus(angle);
            dH.assignMinus(diff.invert().times(dsP));
            dK.assignMinus(diff.invert().times(dsP));
            dH_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            dK_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignDivide(diff);
            K.assignDivide(diff);
            H_corr.assignDivide(sP0.minus(angle));
            K_corr.assignDivide(sP0.minus(angle));
        }
          for (Complex angle : schottky.getAngles()[1]) {
            Complex diff = sP.minus(angle);
            dG.assignMinus(diff.invert().times(dsP));
            dK.assignPlus(diff.invert().times(dsP));
            dG_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            dK_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            G.assignDivide(diff);
            K.assignTimes(diff);
            G_corr.assignDivide(sP0.minus(angle));
            K_corr.assignTimes(sP0.minus(angle));
        }
          for (Complex angle : schottky.getAngles()[2]) {
            Complex diff = sP.minus(angle);
            dH.assignPlus(diff.invert().times(dsP));
            dK.assignMinus(diff.invert().times(dsP));
            dH_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            dK_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignTimes(diff);
            K.assignDivide(diff);
            H_corr.assignTimes(sP0.minus(angle));
            K_corr.assignDivide(sP0.minus(angle));
        }
          for (Complex angle : schottky.getAngles()[3]) {
            Complex diff = sP.minus(angle);
            dG.assignPlus(diff.invert().times(dsP));
            dK.assignPlus(diff.invert().times(dsP));
            dG_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            dK_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            G.assignTimes(diff);
            K.assignTimes(diff);
            G_corr.assignTimes(sP0.minus(angle));
            K_corr.assignTimes(sP0.minus(angle));
        }
    }

    // public Complex aztecArcticCurve(final Complex P, Complex circleCenter, final double accuracy) {
    //   Complex[] diffsP = getDifferentialsQuad(P, accuracy);
    //   Complex[] diffsPDer = {dXi1Der, dXi2Der, dXiAztecDer};
    //   Complex R1 = diffsP[0].divide(diffsP[2]);
    //   Complex R2 = diffsP[1].divide(diffsP[2]);
    //   Complex R3 = diffsP[0].divide(diffsP[1]);
    //   Complex eta = diffsPDer[2].times(diffsP[1]).minus(diffsPDer[1].times(diffsP[2]));
    //   eta.assignDivide(diffsPDer[0].times(diffsP[1]).minus(diffsPDer[1].times(diffsP[0])));
    //   Complex psi = diffsPDer[2].times(diffsP[0]).minus(diffsPDer[0].times(diffsP[2]));
    //   psi.assignDivide(diffsPDer[1].times(diffsP[0]).minus(diffsPDer[0].times(diffsP[1])));
    //   // if (Math.abs(psi.im) > 0.1 | Math.abs(eta.im) > 0.1) {
    //   //   System.out.println("non-real");
    //   // }
    //   return new Complex(psi.re, eta.re);
    // }

    // public Complex aztecArcticCurveReal(final Complex P, final double accuracy) {
    //   // numeric approximation of derivatives approach
    //   Complex[] diffsP = getDifferentialsQuad(P, accuracy);
    //   double epsilon = 0.0001;
    //   Complex PDelta = P.plus(new Complex(epsilon, 0));
    //   Complex[] diffsPDelta = getDifferentialsQuad(PDelta, accuracy);
    //   Complex R1Inv = diffsP[2].divide(diffsP[0]);
    //   Complex R2Inv = diffsP[2].divide(diffsP[1]);
    //   Complex xi = diffsP[0].divide(diffsP[1]);
    //   Complex R1InvDelta = diffsPDelta[2].divide(diffsPDelta[0]);
    //   Complex R2InvDelta = diffsPDelta[2].divide(diffsPDelta[1]);
    //   Complex xiDelta = diffsPDelta[0].divide(diffsPDelta[1]);
    //   double eta = (R2Inv.re - R2InvDelta.re) / (xi.re - xiDelta.re);
    //   double psi = (R1Inv.re - R1InvDelta.re) / (xi.invert().re - xiDelta.invert().re);
    //   return new Complex(psi, eta);
    // }
}
