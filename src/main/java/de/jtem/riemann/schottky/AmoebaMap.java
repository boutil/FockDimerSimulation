package de.jtem.riemann.schottky;
import de.jtem.mfc.field.Complex;

import java.io.Serializable;
import java.lang.Math;

import org.jzy3d.plot3d.pipelines.NotImplementedException;


public class AmoebaMap implements Serializable{

    final SchottkyDimers schottky;

    final int numGenerators;

    int updateID;

    double theta1, q1;

    // measures errors for the elements in some way?
    final double[][] rho;


    AmoebaMap(SchottkyDimers schottky, Complex P0) {
        this.schottky = schottky;
        numGenerators = schottky.numGenerators;
        updateID = schottky.updateID;
        rho = new double[2][numGenerators];
        this.P0.assign(P0);
    }
    
    final Complex zOfRho = new Complex(Double.NaN);
  
    final Complex P = new Complex();
    final Complex P0 = new Complex();

    final Complex A = new Complex();
    final Complex B = new Complex();

    double product_re = 1.0;
    double product_im = 1.0;
    // H for dXi1
    // G for dXi2
    // K for the boundary value differential.
    final Complex H = new Complex();
    final Complex G = new Complex();
    final Complex K = new Complex();

    final Complex dH = new Complex();
    final Complex dG = new Complex();
    final Complex dK = new Complex();

    final Complex dH_Der = new Complex();
    final Complex dG_Der = new Complex();
    final Complex dK_Der = new Complex();

    final Complex H_corr = new Complex();
    final Complex G_corr = new Complex();
    final Complex K_corr = new Complex();
  
    final Complex sP = new Complex();
    final Complex sP0 = new Complex();
    final Complex dsP = new Complex();

    final Complex xi1 = new Complex();
    final Complex xi2 = new Complex();
    final Complex xiBoundary = new Complex();

    final Complex dXi1 = new Complex();
    final Complex dXi2 = new Complex();
    final Complex dXiBoundary = new Complex();

    final Complex dXi1Der = new Complex();
    final Complex dXi2Der = new Complex();
    final Complex dXiBoundaryDer = new Complex();
  
    // final Complex a = new Complex();
    // final Complex b = new Complex();
  
    final Complex d = new Complex();
  
    int n;
  
    long [] noe;
  
    double acc;
    double eps; // = acc / maxNumberOfElements;
    double maxError = 100000;// Double.MAX_VALUE; // set this dynmically in a sensible way.
    
    void update() {
  
      if (updateID == schottky.updateID) {
        return;
      }
  
      updateID = schottky.updateID;
  
      zOfRho.assign(Double.NaN);
  
      theta1 = schottky.theta1;
      q1 = schottky.q1;
      }

  /**
   * Returns rhoIntegral( sigma, z ) for zOfRho.
   * Thus you have to call perepareRho( z ) first.
   * @param sigma
   * @return rho( sigma, zOfRhoDiff )
   */
  final double rho(SchottkyGroupElement sigma) {
    return rho[sigma.leftIsInvert][sigma.left];
  }

  /**
   * Prepares rho for value z.
   * @param z
   */
  final void prepareRho(Complex z) {

    if (z.equals(zOfRho)) {
      return;
    }

    zOfRho.assign(z);

    if (schottky.useFancyError) {
      final double q_ = theta1 * theta1;
      final double R2 = schottky.r(q_);
      final double R2Minus = schottky.rMinus(R2, q_);
      final double R2Plus = schottky.rPlus(R2, q_);

      //final double[][] k2 = schottky.k2(z);
      final double[][] k2 = schottky.kIndexed(z, 3);

      for (int j = 0; j < 2; j++) {
        for (int m = 0; m < numGenerators; m++) {

          double rhoJM = 0;

          for (int i = 0; i < 2; i++) {
            for (int n = 0; n < numGenerators; n++) {

              if (n == m) {
                if (i == j) {
                  rhoJM += R2Plus / k2[i][n];
                }
                else {
                  rhoJM += R2Minus / k2[i][n];
                }
              }
              else {
                rhoJM += R2 / k2[i][n];
              }
            }
          }
          rho[j][m] = rhoJM;
        }
      }
    }
    else {


      //double dOfZ = schottky.d2(z);
      double dOfZ = schottky.k(z,3);

      for (int j = 0; j < 2; j++) {
        for (int m = 0; m < numGenerators; m++) {
          rho[j][m] = 1 / dOfZ
              / (1 - q1);
        }
      }
    }
  }

  protected void cleanIncrements() {
      dH.assign(0);
      dG.assign(0);
      dK.assign(0);
      dH_Der.assign(0);
      dG_Der.assign(0);
      dK_Der.assign(0);
      H.assign(1);
      G.assign(1);
      K.assign(1);
      H_corr.assign(1);
      G_corr.assign(1);
      K_corr.assign(1);
  }

  protected void cleanDiffs() {
    xi1.assign(0, 0);
    xi2.assign(0, 0);
    xiBoundary.assign(0, 0);

    dXi1.assign(0, 0);
    dXi2.assign(0, 0);
    dXiBoundary.assign(0, 0);

    dXi1Der.assign(0, 0);
    dXi2Der.assign(0, 0);
    dXiBoundaryDer.assign(0, 0);
  }

      // // Adds the U1 correction to the Amoeba map. g=1 only for now. Therefore n = 0 only instead of doing some Matrix inversion.

      // final void addU1Correction(Complex P) {
      //   // double abelIntegral = schottky.abelianIntegralOf1stKind(P, 0).re;
      //   Complex P0 = new Complex(0.5, 0);
      //   double abelIntegral = Math.log(Complex.crossRatio(B, P0, A, P).abs());
      //   Complex factor1 = new Complex();
      //   Complex factor2 = new Complex();
      //   factor1.assignCrossRatio(schottky.getA(0), alpha, schottky.getB(0), gamma);
      //   double f1 = Math.log(factor1.re) / Math.log(schottky.getMu(0).re) * abelIntegral;
      //   factor2.assignCrossRatio(schottky.getA(0), beta, schottky.getB(0), gamma);
      //   double f2 = Math.log(factor2.re) / Math.log(schottky.getMu(0).re) * abelIntegral;
      //   product_re -= f1;
      //   product_im -= f2;
      // }



    // calculates the amoebaMap based on a quadGrid setup. Thus there need to be four angles picked.
    final void processSchottkyElement(final SchottkyGroupElement element) {

      if (element.updateID != updateID) {
        schottky.updateElement(element);
      }
    // Not sure what to do here instead...
      element.diff(B, A, d);
      double error = maxError - 1;
      if (element != schottky.id) {
        error = (Math.abs(d.re) + Math.abs(d.im)) * rho(element);
      }

      if (error > maxError || Double.isNaN(error)) {
        return;
      }
  
      // This way we get genus 0.
      // if (element.wordLength > 0) {
      //   return;
      // }

      if (error * noe[element.wordLength] < acc || error < eps ) {
        acc += acc / noe[element.wordLength] - error;
        return;
      }

      cleanIncrements();

      
      element.applyTo(P, sP);
      element.applyTo(P0, sP0);

      element.applyDifferentialTo(P, dsP);

      calculateIncrements(element);

      addIncrements();

      if (element.child == null) {
        schottky.createLeftChilds(element);
      }
  
      final SchottkyGroupElement[] child = element.child;
  
      final int numChildren = child.length;
  
      String[] childrenWords = new String[numChildren];
      for (int i = 0; i < numChildren; i++) {
        childrenWords[i] = child[i].word();
      }
      for (int i = 0; i < numChildren; i++) {
        // String childWord = child[i].word();
        processSchottkyElement(child[i]);
      }
    }

    protected void calculateIncrements(final SchottkyGroupElement element) {
      throw new NotImplementedException();
    }

    protected void addIncrements() {
      if(!dH.isNaN()){
        dXi1.assignPlus(dH);
      }
      if(!dG.isNaN()){
        dXi2.assignPlus(dG);
      }
      if(!dK.isNaN()){
        dXiBoundary.assignPlus(dK);
      }
      if(!H.divide(H_corr).log().isNaN()){
        xi1.assignPlus(H.divide(H_corr).log());
      }
      if(!G.divide(G_corr).log().isNaN()){
        xi2.assignPlus(G.divide(G_corr).log());
      }
      if(!K.divide(K_corr).log().isNaN()){
        xiBoundary.assignPlus(K.divide(K_corr).log());
      }
      if(!dH_Der.isNaN()){
        dXi1Der.assignPlus(dH_Der);
      }
      if(!dG_Der.isNaN()){
        dXi2Der.assignPlus(dG_Der);
      }
      if(!dK_Der.isNaN()){
        dXiBoundaryDer.assignPlus(dK_Der);
      }
    }

      final Complex[] getDifferentials(final Complex P, final double accuracy) {
    
        // calculates and returns xi1, xi2 and xiTilde.
        // (xi1.re, xi2.re) is then Amoeba Map.
        this.noe = schottky.numOfElementsWithWordLength;
        
        this.A.assign(schottky.fixpoint[n][0]);
        this.B.assign(schottky.fixpoint[n][1]);

        this.P.assign(P);

        cleanDiffs();

        this.acc = accuracy;
        this.eps = accuracy / schottky.maxNumOfElements;

        prepareRho(P);

        processSchottkyElement(schottky.id);

        // Correction in the U1 uniformization case.
        // if(schottky.uniformization == 1) {
        //   addU1Correction(P);
        // }

        if(this.acc < 0 ) // this test is needed because of the eps crieteria
          throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
        
        // r.assign(new Complex(Math.log(product_re), Math.log(product_im)));
        return new Complex[]{dXi1.copy(), dXi2.copy(), dXiBoundary.copy(), xi1.copy(), xi2.copy(), xiBoundary.copy(), dXi1Der.copy(), dXi2Der.copy(), dXiBoundaryDer.copy()};
      }

      public Complex amoebaMap(final Complex P, final double accuracy) {
        Complex[] diffs = getDifferentials(P, accuracy);
        return new Complex(diffs[3].re, diffs[4].re);
      }
      
      
      public Complex boundaryMap(final Complex P, final double accuracy) {
        Complex[] diffs = getDifferentials(P, accuracy);

        // Complex R1 = diffs[0].divide(diffs[2]);
        // Complex R2 = diffs[1].divide(diffs[2]);
        // psi, eta are the aztec diamond coordinates.
        // if ((Math.abs(R1.im/R2.im) > 50 || Math.abs(R2.im/R1.im) > 50) && P.im > 0.01) {
        //   System.out.println(R1.im/R2.im);
        // }
        // double psi = - R2.re + R1.re * (R2.im/R1.im);
        // double eta = R1.re - R2.re * (R1.im/R2.im);
        // double psi = -R2.invert().im / R1.divide(R2).im;
        // double eta = R1.invert().im / R2.divide(R1).im;
        // return new Complex(psi, -eta);

        // probably should switch to a more numerically stable version here. Consider using BigDecimal.
        double factor = diffs[0].times(diffs[1].conjugate()).im;
        double psi = diffs[0].times(diffs[2].conjugate()).im;
        double eta = diffs[1].times(diffs[2].conjugate()).im;

        // if (Math.abs(psi/factor - eta/factor) > 1) {
        //   System.out.println(psi/factor + ", " + eta/factor);
        // }
        return new Complex(psi/factor, -eta/factor);
      }

}
