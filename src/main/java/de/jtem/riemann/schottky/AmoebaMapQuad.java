package de.jtem.riemann.schottky;

import de.jtem.mfc.field.Complex;

public class AmoebaMapQuad extends AmoebaMap{

    AmoebaMapQuad(SchottkyDimers schottky, Complex P0) {
        super(schottky, P0);
    }

    AmoebaMapQuad(SchottkyDimers schottky, Complex P0, double[] boundaryResidues) {
      this(schottky, P0);
      this.boundaryResidues = boundaryResidues;
  }

    @Override
    protected void calculateIncrements(final SchottkyGroupElement element) {
        for (Complex angle : schottky.getAngles()[0]) {
            Complex diff = P.minus(element.applyTo(angle));
            dH.assignMinus(diff.invert());
            dK.assignMinus(diff.invert());
            dH_Der.assignPlus(diff.pow(2).invert());
            dK_Der.assignPlus(diff.pow(2).invert());
            H.assignMinus(diff.log());
            K.assignMinus(diff.log());
            H_corr.assignPlus(P0.minus(element.applyTo(angle)).log());
            K_corr.assignPlus(P0.minus(element.applyTo(angle)).log());
        }
          for (Complex angle : schottky.getAngles()[1]) {
            Complex diff = P.minus(element.applyTo(angle));
            dG.assignMinus(diff.invert());
            dK.assignPlus(diff.invert());
            dG_Der.assignPlus(diff.pow(2).invert());
            dK_Der.assignMinus(diff.pow(2).invert());
            G.assignMinus(diff.log());
            K.assignPlus(diff.log());
            G_corr.assignMinus(P0.minus(element.applyTo(angle)).log());
            K_corr.assignPlus(P0.minus(element.applyTo(angle)).log());
        }
          for (Complex angle : schottky.getAngles()[2]) {
            Complex diff = P.minus(element.applyTo(angle));
            dH.assignPlus(diff.invert());
            dK.assignMinus(diff.invert());
            dH_Der.assignMinus(diff.pow(2).invert());
            dK_Der.assignPlus(diff.pow(2).invert());
            H.assignPlus(diff.log());
            K.assignMinus(diff.log());
            H_corr.assignMinus(P0.minus(element.applyTo(angle)).log());
            K_corr.assignPlus(P0.minus(element.applyTo(angle)).log());
        }
          for (Complex angle : schottky.getAngles()[3]) {
            Complex diff = P.minus(element.applyTo(angle));
            dG.assignPlus(diff.invert());
            dK.assignPlus(diff.invert());
            dG_Der.assignMinus(diff.pow(2).invert());
            dK_Der.assignMinus(diff.pow(2).invert());
            G.assignPlus(diff.log());
            K.assignPlus(diff.log());
            G_corr.assignMinus(P0.minus(element.applyTo(angle)).log());
            K_corr.assignMinus(P0.minus(element.applyTo(angle)).log());
        }
    }



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
