package de.jtem.riemann.schottky;
import de.jtem.mfc.field.Complex;
import java.lang.Math;
import java.util.Objects;


public class AmoebaMap {

    final SchottkyDimers schottky;

    final int numGenerators;

    int updateID;

    double theta1, q1;

    // measures errors for the elements in some way?
    final double[][] rho;

    Complex alpha;
    Complex beta;
    Complex gamma;


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

    final Complex H_corr = new Complex();
    final Complex G_corr = new Complex();
  
    final Complex sP = new Complex();
    final Complex sP0 = new Complex();
  
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



}
