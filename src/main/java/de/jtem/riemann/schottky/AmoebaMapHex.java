package de.jtem.riemann.schottky;

import de.jtem.mfc.field.Complex;

public class AmoebaMapHex extends AmoebaMap{
    
    boolean isDoubleCover;

    AmoebaMapHex(SchottkyDimers schottky, Complex P0) {
        super(schottky, P0);

        isDoubleCover = (schottky.getAngles().length == 6);
    }


    @Override
    protected void calculateIncrements(final SchottkyGroupElement element) {
        if (isDoubleCover) {
            doubleCoverIncrements(element);
        } else {
            simpleIncrements(element);
        }
    }

    private void simpleIncrements(final SchottkyGroupElement element) {
        for (Complex angle : schottky.getAngles()[0]) {
            Complex diff = sP.minus(angle);
            dH.assignPlus(diff.invert().times(dsP));
            dH_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignTimes(diff);
            H_corr.assignTimes(sP0.minus(angle));
        }
          for (Complex angle : schottky.getAngles()[1]) {
            Complex diff = sP.minus(angle);
            dG.assignPlus(diff.invert().times(dsP));
            dG_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            G.assignTimes(diff);
            G_corr.assignTimes(sP0.minus(angle));
        }
          for (Complex angle : schottky.getAngles()[2]) {
            Complex diff = sP.minus(angle);
            dH.assignMinus(diff.invert().times(dsP));
            dG.assignMinus(diff.invert().times(dsP));
            dH_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            dG_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignDivide(diff);
            G.assignDivide(diff);
            H_corr.assignDivide(sP0.minus(angle));
            G_corr.assignDivide(sP0.minus(angle));
        }
    }

    private void doubleCoverIncrements(final SchottkyGroupElement element) {
        for (Complex angle : schottky.getAngles()[0]) {
            Complex diff = sP.minus(angle);
            dH.assignPlus(diff.invert().times(dsP));
            dH_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignTimes(diff);
            H_corr.assignTimes(sP0.minus(angle));
            dK.assignPlus(diff.invert().times(dsP));
        }
          for (Complex angle : schottky.getAngles()[1]) {
            Complex diff = sP.minus(angle);
            dG.assignPlus(diff.invert().times(dsP));
            dG_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            G.assignTimes(diff);
            G_corr.assignTimes(sP0.minus(angle));
            dK.assignMinus(diff.invert().times(dsP));
        }
        for (Complex angle : schottky.getAngles()[2]) {
            Complex diff = sP.minus(angle);
            dH.assignMinus(diff.invert().times(dsP));
            dG.assignMinus(diff.invert().times(dsP));
            dH_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            dG_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignDivide(diff);
            G.assignDivide(diff);
            H_corr.assignDivide(sP0.minus(angle));
            G_corr.assignDivide(sP0.minus(angle));
            dK.assignPlus(diff.invert().times(dsP));
        }
        for (Complex angle : schottky.getAngles()[3]) {
            Complex diff = sP.minus(angle);
            dH.assignPlus(diff.invert().times(dsP));
            dH_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignTimes(diff);
            H_corr.assignTimes(sP0.minus(angle));
            dK.assignMinus(diff.invert().times(dsP));
        }
        for (Complex angle : schottky.getAngles()[4]) {
            Complex diff = sP.minus(angle);
            dG.assignPlus(diff.invert().times(dsP));
            dG_Der.assignPlus(P.minus(element.applyTo(angle)).pow(2).invert());
            G.assignTimes(diff);
            G_corr.assignTimes(sP0.minus(angle));
            dK.assignPlus(diff.invert().times(dsP));
        }
        for (Complex angle : schottky.getAngles()[5]) {
            Complex diff = sP.minus(angle);
            dH.assignMinus(diff.invert().times(dsP));
            dG.assignMinus(diff.invert().times(dsP));
            dH_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            dG_Der.assignMinus(P.minus(element.applyTo(angle)).pow(2).invert());
            H.assignDivide(diff);
            G.assignDivide(diff);
            H_corr.assignDivide(sP0.minus(angle));
            G_corr.assignDivide(sP0.minus(angle));
            dK.assignMinus(diff.invert().times(dsP));
        }
    }
}
