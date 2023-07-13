package de.jtem.riemann.schottky;
import de.jtem.mfc.field.Complex;

import java.io.Serializable;
import java.lang.Math;
import java.util.Objects;


public class AmoebaMap implements Serializable{

    final SchottkyDimers schottky;

    final int numGenerators;

    int updateID;

    double theta1, q1;

    // measures errors for the elements in some way?
    final double[][] rho;

    // For hex grid:
    Complex alpha;
    Complex beta;
    Complex gamma;

    // For quad lattice assumed to be alpha-, beta-, alpha+, beta+.
    // TODO For hex grid assumed to be alpha, beta, gamma. 
    Complex[] angles;


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
    final Complex H = new Complex();
    final Complex G = new Complex();
    final Complex K = new Complex();

    final Complex dH = new Complex();
    final Complex dG = new Complex();
    final Complex dK = new Complex();

    final Complex H_corr = new Complex();
    final Complex G_corr = new Complex();
    final Complex K_corr = new Complex();
  
    final Complex sP = new Complex();
    final Complex sP0 = new Complex();
    final Complex dsP = new Complex();

    final Complex xi1 = new Complex();
    final Complex xi2 = new Complex();
    final Complex xiAztec = new Complex();

    final Complex dXi1 = new Complex();
    final Complex dXi2 = new Complex();
    final Complex dXiAztec = new Complex();

    final Complex dXi1Der = new Complex();
    final Complex dXi2Der = new Complex();
    final Complex dXiAztecDer = new Complex();
  
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
    // calculates the amoebaMap based on a hexGrid setup. Thus there need to be three angles picked.
    final void hexGrid(final SchottkyGroupElement element) {

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
    
        if (element.wordLength > 0) {
          return;
        }

        if (error * noe[element.wordLength] < acc || error < eps ) {
          acc += acc / noe[element.wordLength] - error;
          return;
        }

        
        element.applyTo(P, sP);
        element.applyTo(P0, sP0);
        
        H.assignDivide(sP.minus(alpha), sP.minus(gamma));
        G.assignDivide(sP.minus(beta), sP.minus(gamma));
        H_corr.assignDivide(sP0.minus(gamma), sP0.minus(alpha));
        G_corr.assignDivide(sP0.minus(gamma), sP0.minus(beta));

        product_re += Math.log(H.abs() * H_corr.abs());
        product_im += Math.log(G.abs() * G_corr.abs());

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
          hexGrid(child[i]);
        }
      }
    
      final void hexGrid
          (final Complex r,
           final Complex P,
           final Complex alpha, final Complex beta, final Complex gamma,
           final double accuracy) {
    
        this.noe = schottky.numOfElementsWithWordLength;
        
        this.A.assign(schottky.fixpoint[n][0]);
        this.B.assign(schottky.fixpoint[n][1]);

        this.P.assign(P);
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;

        this.acc = accuracy;
        this.eps = accuracy / schottky.maxNumOfElements;

        prepareRho(P);
        // H.assignDivide(P.minus(alpha), P.minus(gamma));
        // G.assignDivide(P.minus(beta), P.minus(gamma));
        // product_re = H.abs();
        // product_im = G.abs();
        product_re = 0;
        product_im = 0;

        hexGrid(schottky.id);

        // Correction in the U1 uniformization case.
        if(schottky.uniformization == 1) {
          addU1Correction(P);
        }
        // for (int i = 0; i < numGenerators; i++) {
        //   hexGrid(schottky.generator[i]);
        //   hexGrid(schottky.generatorInv[i]);
        // }

        if(this.acc < 0 ) // this test is needed because of the eps crieteria
          throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
        
        // r.assign(new Complex(Math.log(product_re), Math.log(product_im)));
        r.assign(new Complex(product_re, product_im));
      }

      // Adds the U1 correction to the Amoeba map. g=1 only for now. Therefore n = 0 only instead of doing some Matrix inversion.
      final void addU1Correction(Complex P) {
        // double abelIntegral = schottky.abelianIntegralOf1stKind(P, 0).re;
        Complex P0 = new Complex(0.5, 0);
        double abelIntegral = Math.log(Complex.crossRatio(B, P0, A, P).abs());
        Complex factor1 = new Complex();
        Complex factor2 = new Complex();
        factor1.assignCrossRatio(schottky.getA(0), alpha, schottky.getB(0), gamma);
        double f1 = Math.log(factor1.re) / Math.log(schottky.getMu(0).re) * abelIntegral;
        factor2.assignCrossRatio(schottky.getA(0), beta, schottky.getB(0), gamma);
        double f2 = Math.log(factor2.re) / Math.log(schottky.getMu(0).re) * abelIntegral;
        product_re -= f1;
        product_im -= f2;
      }



    // calculates the amoebaMap based on a quadGrid setup. Thus there need to be four angles picked.
    final void quadGrid(final SchottkyGroupElement element) {

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
  
      // if (element.wordLength > 0) {
      //   return;
      // }

      if (error * noe[element.wordLength] < acc || error < eps ) {
        acc += acc / noe[element.wordLength] - error;
        return;
      }

      
      element.applyTo(P, sP);
      element.applyTo(P0, sP0);

      element.applyDifferentialTo(P, dsP);


      dH.assign(sP.minus(angles[2]).invert().minus(sP.minus(angles[0]).invert()).times(dsP));
      dG.assign(sP.minus(angles[3]).invert().minus(sP.minus(angles[1]).invert()).times(dsP));
      dK.assign(sP.minus(angles[1]).invert().plus(sP.minus(angles[3]).invert()).minus(sP.minus(angles[2]).invert()).minus(sP.minus(angles[0]).invert()).times(dsP));



      H.assignDivide(sP.minus(angles[2]), sP.minus(angles[0]));
      G.assignDivide(sP.minus(angles[3]), sP.minus(angles[1]));
      K.assignDivide(sP.minus(angles[1]).times(sP.minus(angles[3])), sP.minus(angles[0]).times(sP.minus(angles[2])));
      // K.assignDivide(sP.minus(angles[2]).times(sP.minus(angles[0])), sP.minus(angles[1]).times(sP.minus(angles[3])));
      H_corr.assignDivide(sP0.minus(angles[2]), sP0.minus(angles[0]));
      G_corr.assignDivide(sP0.minus(angles[3]), sP0.minus(angles[1]));
      K_corr.assignDivide(sP0.minus(angles[1]).times(sP0.minus(angles[3])), sP0.minus(angles[2]).times(sP0.minus(angles[0])));
      // K_corr.assignDivide(sP0.minus(angles[2]).times(sP0.minus(angles[0])), sP0.minus(angles[1]).times(sP0.minus(angles[3])));

      if(!dH.isNaN()){
        dXi1.assignPlus(dH);
      }
      if(!dG.isNaN()){
        dXi2.assignPlus(dG);
      }
      if(!dK.isNaN()){
        dXiAztec.assignPlus(dK);
      }
      if(!H.divide(H_corr).log().isNaN()){
        xi1.assignPlus(H.divide(H_corr).log());
      }
      if(!G.divide(G_corr).log().isNaN()){
        xi2.assignPlus(G.divide(G_corr).log());
      }
      if(!K.divide(K_corr).log().isNaN()){
        xiAztec.assignPlus(K.divide(K_corr).log());
      }

      // add boundary data correction
      // K.assignDivide(sP.minus(angles[]))

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
        quadGrid(child[i]);
      }
    }

    final void quadGridCircleDers(final SchottkyGroupElement element, Complex circleCenter) {

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
  
      // if (element.wordLength > 0) {
      //   return;
      // }

      if (error * noe[element.wordLength] < acc || error < eps ) {
        acc += acc / noe[element.wordLength] - error;
        return;
      }

      dH.assign(P.minus(element.applyTo(angles[2])).pow(2).invert().minus(P.minus(element.applyTo(angles[0])).pow(2).invert()));
      dG.assign(P.minus(element.applyTo(angles[3])).pow(2).invert().minus(P.minus(element.applyTo(angles[1])).pow(2).invert()));
      dK.assign(P.minus(element.applyTo(angles[1])).pow(2).invert().plus(sP.minus(element.applyTo(angles[3])).pow(2).invert()).minus(sP.minus(element.applyTo(angles[2])).pow(2).invert()).minus(sP.minus(element.applyTo(angles[0])).pow(2).invert()));

      if(!dH.isNaN()){
        dXi1Der.assignPlus(dH);
      }
      if(!dG.isNaN()){
        dXi2Der.assignPlus(dG);
      }
      if(!dK.isNaN()){
        dXiAztecDer.assignPlus(dK);
      }

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
        quadGridCircleDers(child[i], circleCenter);
      }
    }

      final Complex[] getDifferentials
          (final Complex P,
           final Complex[] angles,
           final double accuracy) {
    
        // calculates and returns xi1, xi2 and xiTilde.
        // (xi1.re, xi2.re) is then Amoeba Map.
        this.noe = schottky.numOfElementsWithWordLength;
        
        this.A.assign(schottky.fixpoint[n][0]);
        this.B.assign(schottky.fixpoint[n][1]);

        this.P.assign(P);
        this.angles = angles;

        xi1.assign(0, 0);
        xi2.assign(0, 0);
        xiAztec.assign(0, 0);

        dXi1.assign(0, 0);
        dXi2.assign(0, 0);
        dXiAztec.assign(0, 0);

        this.acc = accuracy;
        this.eps = accuracy / schottky.maxNumOfElements;

        prepareRho(P);

        quadGrid(schottky.id);

        // Correction in the U1 uniformization case.
        if(schottky.uniformization == 1) {
          addU1Correction(P);
        }

        if(this.acc < 0 ) // this test is needed because of the eps crieteria
          throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
        
        // r.assign(new Complex(Math.log(product_re), Math.log(product_im)));
        return new Complex[]{dXi1.copy(), dXi2.copy(), dXiAztec.copy(), xi1.copy(), xi2.copy(), xiAztec.copy()};
      }

      public Complex[] getDifferentialDers(final Complex P, final Complex angles[], Complex circleCenter, final double accuracy) {
        this.noe = schottky.numOfElementsWithWordLength;
        
        this.A.assign(schottky.fixpoint[n][0]);
        this.B.assign(schottky.fixpoint[n][1]);

        this.P.assign(P);
        this.angles = angles;

        dXi1Der.assign(0, 0);
        dXi2Der.assign(0, 0);
        dXiAztecDer.assign(0, 0);

        
        this.acc = accuracy;
        this.eps = accuracy / schottky.maxNumOfElements;
        
        prepareRho(P);
        
        quadGridCircleDers(schottky.id, circleCenter);

        dXi1Der.assignTimes(P.minus(circleCenter).timesI().neg());
        dXi2Der.assignTimes(P.minus(circleCenter).timesI().neg());
        dXiAztecDer.assignTimes(P.minus(circleCenter).timesI().neg());

        if(this.acc < 0 ) // this test is needed because of the eps crieteria
          throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
        
        // r.assign(new Complex(Math.log(product_re), Math.log(product_im)));
        return new Complex[]{dXi1Der.copy(), dXi2Der.copy(), dXiAztecDer.copy()};
      }

      public Complex amoebaMapQuadGrid(final Complex P, final Complex angles[], final double accuracy) {
        Complex[] diffs = getDifferentials(P, angles, accuracy);
        return new Complex(diffs[3].re, diffs[4].re);
      }
      
      public Complex aztecMap(final Complex P, final Complex angles[], final double accuracy) {
        Complex[] diffs = getDifferentials(P, angles, accuracy);
        // Complex R1 = diffs[0].divide(diffs[2]);
        // Complex R2 = diffs[1].divide(diffs[2]);
        // psi, eta are the aztec diamond coordinates.
        // if ((Math.abs(R1.im/R2.im) > 50 || Math.abs(R2.im/R1.im) > 50) && P.im > 0.01) {
        //   System.out.println(R1.im/R2.im);
        // }
        // double psi = - R2.re + R1.re * (R2.im/R1.im);
        // double eta = R1.re - R2.re * (R1.im/R2.im);
        // return new Complex(1/psi, 1/eta);
        // double psi = -R2.invert().im / R1.divide(R2).im;
        // double eta = R1.invert().im / R2.divide(R1).im;
          double factor = diffs[0].times(diffs[1].conjugate()).im;
          double psi = diffs[0].times(diffs[2].conjugate()).im;
          double eta = diffs[1].times(diffs[2].conjugate()).im;
        return new Complex(psi/factor, -eta/factor);
      }

      public Complex aztecArcticCurve(final Complex P, final Complex[] angles, Complex circleCenter, final double accuracy) {
        Complex[] diffsP = getDifferentials(P, angles, accuracy);
        Complex[] diffsPDer = getDifferentialDers(P, angles, circleCenter, accuracy);
        Complex R1 = diffsP[0].divide(diffsP[2]);
        Complex R2 = diffsP[1].divide(diffsP[2]);
        Complex R3 = diffsP[0].divide(diffsP[1]);
        Complex eta = diffsPDer[2].times(diffsP[1]).minus(diffsPDer[1].times(diffsP[2]));
        eta.assignDivide(diffsPDer[0].times(diffsP[1]).minus(diffsPDer[1].times(diffsP[0])));
        Complex psi = diffsPDer[2].times(diffsP[0]).minus(diffsPDer[0].times(diffsP[2]));
        psi.assignDivide(diffsPDer[1].times(diffsP[0]).minus(diffsPDer[0].times(diffsP[1])));
        // if (Math.abs(psi.im) > 0.1 | Math.abs(eta.im) > 0.1) {
        //   System.out.println("non-real");
        // }
        return new Complex(psi.re, eta.re);
      }

      public Complex aztecArcticCurveReal(final Complex P, final Complex[] angles, final double accuracy) {
        Complex[] diffsP = getDifferentials(P, angles, accuracy);
        double epsilon = 0.0001;
        Complex PDelta = P.plus(new Complex(epsilon, 0));
        Complex[] diffsPDelta = getDifferentials(PDelta, angles, accuracy);
        Complex R1Inv = diffsP[2].divide(diffsP[0]);
        Complex R2Inv = diffsP[2].divide(diffsP[1]);
        Complex xi = diffsP[0].divide(diffsP[1]);
        Complex R1InvDelta = diffsPDelta[2].divide(diffsPDelta[0]);
        Complex R2InvDelta = diffsPDelta[2].divide(diffsPDelta[1]);
        Complex xiDelta = diffsPDelta[0].divide(diffsPDelta[1]);
        double eta = (R2Inv.re - R2InvDelta.re) / (xi.re - xiDelta.re);
        double psi = (R1Inv.re - R1InvDelta.re) / (xi.invert().re - xiDelta.invert().re);
        return new Complex(psi, eta);
      }

}
