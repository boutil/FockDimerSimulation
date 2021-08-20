/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2009, Technische Universität Berlin, jTEM
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

-	Redistributions of source code must retain the above copyright notice, 
	this list of conditions and the following disclaimer.

-	Redistributions in binary form must reproduce the above copyright notice, 
	this list of conditions and the following disclaimer in the documentation 
	and/or other materials provided with the distribution.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
OF SUCH DAMAGE.
**/

package de.jtem.riemann.schottky;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * @author Markus Schmies
 * @version 1.0
 */
class AbelianDifferential
    implements Serializable {

  final Schottky schottky;

  final double [][] rho;

  final double[] L1;

  final int numGenerators;

  int updateID;

  double kappa2, q2, f2, d2ForZero, maxInIsometricCircles;

  static final Complex ZERO = new Complex();

  AbelianDifferential(Schottky schottky) {
    this.schottky = schottky;
    numGenerators = schottky.numGenerators;
    updateID = schottky.updateID;
    rho = new double[2][numGenerators];
    L1 = new double[numGenerators];
  }

  final Complex zOfRho
      = new Complex(Double.NaN);

  final Complex z = new Complex();
  final Complex s = new Complex();
  final Complex d = new Complex();
  final Complex H = new Complex();

  final Complex A = new Complex();
  final Complex B = new Complex();

  final Complex a = new Complex();
  final Complex b = new Complex();

  int n;

  long [] noe;

  double acc;
  double eps; // = acc / maxNumberOfElements;
  
  void update() {

    if (updateID == schottky.updateID)
      return;

    updateID = schottky.updateID;

    if (numGenerators > 1) {
      kappa2 = kappa2();

      if( kappa2<0)
        throw new RuntimeException( "kappa is negative?????" );

      q2 = (2 * numGenerators - 1) / kappa2 / kappa2;

      if (q2 >= 1 ) {

        kappa2 = kappa(3);

        if( kappa2<0)
          throw new RuntimeException( "kappa is negative?????" );

        q2 = (2 * numGenerators - 1) / kappa2 / kappa2;


        if (q2 >= 1 ) {
          throw new RuntimeException( "cannot guarantee convergence of poincaree theta series" );
        }
      }



      f2 = 1 / (1-q2);

    }

    zOfRho.assign(Double.NaN);

    for (int n = 0; n < numGenerators; n++) {
      L1[n] = L1(n);
    }
  }

  final double getKappa2() {
    update();
    return kappa2;
  }

  final double getQ2() {
    update();
    return q2;
  }

  /**
   * Returns rho( sigma, z ) for zOfRhoDiff.
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
      final double q_ = 1 / kappa2 / kappa2;
      final double R2 = schottky.r(q_);
      final double R2Minus = schottky.rMinus(R2, q_);
      final double R2Plus = schottky.rPlus(R2, q_);

      //final double[][] k2 = schottky.k2(z);
      final double[][] k2 = schottky.kIndexed(z,3);

      for (int j = 0; j < 2; j++) {
        for (int m = 0; m < numGenerators; m++) {

          double rhoJM = 0;

          for (int i = 0; i < 2; i++) {
            for (int n = 0; n < numGenerators; n++) {

              double d2Sqr = k2[i][n];
              d2Sqr *= d2Sqr;

              if (n == m) {
                if (i == j) {
                  rhoJM += R2Plus / d2Sqr;
                }
                else {
                  rhoJM += R2Minus / d2Sqr;
                }
              }
              else {
                rhoJM += R2 / d2Sqr;
              }
            }
          }
          rho[j][m] = rhoJM;
        }
      }

    }
    else {

      //double dOfZSqr = schottky.d2(z);
      double dOfZSqr = schottky.k(z,3);
      dOfZSqr *= dOfZSqr;

      for (int j = 0; j < 2; j++) {
        for (int m = 0; m < numGenerators; m++) {
          rho[j][m] = 1 / dOfZSqr / (1 - q2);
        }
      }
    }
  }

  final double L(SchottkyGroupElement sigma, double[] LOfASigma) {

    double maxL = 0;

    for (int n = 0; n < numGenerators; n++) {
      if (n != sigma.left) {

        LOfASigma[n] = L(sigma, n);

        if (LOfASigma[n] > maxL) {
          maxL = LOfASigma[n];
        }
      }
      else {
        LOfASigma[n] = 0;
      }
    }

    return maxL;
  }

  final double L(SchottkyGroupElement sigma, int n) {

    if (sigma.left == n) {
      throw new RuntimeException("do not compute L for this n");
    }

    return L(sigma, schottky.fixpoint[n][0], schottky.fixpoint[n][1]);
  }

  final double L(SchottkyGroupElement sigma, Complex A, Complex B) {
    return A.dist(B) / schottky.dist(sigma, A) / schottky.dist(sigma, B);
    }

  final double L1(int n) {

    double max = 0;

    for (int i = 0; i < numGenerators; i++) {

      if (i != n) {

        final double L1OfI = Math.max(L(schottky.generator[i], n),
                                      L(schottky.generatorInv[i], n));

        if (L1OfI > max) {
          max = L1OfI;
        }
      }
    }

    return max;
  }

  final double kappaLBar( SchottkyGroupElement sigma ) {

    SchottkyGroupElement tau = sigma.parent;

    final int i    = sigma.leftIsInvert;
    final int l    = sigma.left;

    final Complex An = schottky.fixpoint[l][0];
    final Complex Bn = schottky.fixpoint[l][1];

    final double sqrtOfAbsMu;

    if (i == 1) {
      sqrtOfAbsMu = 1 / Math.sqrt(schottky.mu[l].abs());
    }
    else {
      sqrtOfAbsMu = Math.sqrt(schottky.mu[l].abs());
    }

    final double v1 = schottky.k(tau, An) / sqrtOfAbsMu
        - schottky.K(tau, Bn) * sqrtOfAbsMu;
    final double v2 = schottky.k(tau, Bn) * sqrtOfAbsMu
        - schottky.K(tau, An) / sqrtOfAbsMu;

    return v1 > v2 ? v1 / An.dist(Bn) : v2 / An.dist(Bn);
  }


      final private double evalKappaLBar(SchottkyGroupElement element,
                                         double kappa, int wordLength ) {

        if (element.wordLength > wordLength )
          return kappa;

        if (element.updateID != updateID) {
          schottky.updateElement(element);
        }

        if( element.wordLength == wordLength ) {
          kappa = Math.min( kappaLBar( element ), kappa );
        }

        if (element.child == null) {
          schottky.createLeftChilds(element);
        }

        final SchottkyGroupElement[] child = element.child;

        final int numOfChilds = child.length;

        for (int i = 0; i < numOfChilds; i++) {
          kappa = evalKappaLBar(child[i], kappa, wordLength );
        }
        return kappa;
      }

      /**
       * Computes kappa for given word lengths.
       * kappa is only defined for word lengths greater then one."
       */
      public final double kappa( int wordLength ) {

        if( wordLength < 2 )
          throw new IllegalArgumentException
              ( "kappa only defined for word legnth greater then one" );

        if (updateID != schottky.updateID)
          update();

        return evalKappaLBar(schottky.id, Double.MAX_VALUE, wordLength );
      }


    /**
     * return kappaRBar( sigma ) with sigma = g(j,m) g(i,n).
     * The computation is not quite right. The definition relates
     * kappaRBar( sigma ) to
     * kappaR( sigma ) which gives the right choice of the two
     * possible values. We just take the maximum of
     * the two possible values, we believe that this has no effect.
     * A possible error would be on the save side and only deliver
     * worse estimates.
     */
    final double kappaRBar( final int j, final int m, final int i, final int r) {

      final int jInv = j==0 ? 1 : 0;

      final Complex An = schottky.fixpoint[r][0];
      final Complex Bn = schottky.fixpoint[r][1];

      final double sqrtOfAbsMu;

      if (i == 1) {
        sqrtOfAbsMu = 1 / Math.sqrt(schottky.mu[r].abs());
      } else {
        sqrtOfAbsMu =     Math.sqrt(schottky.mu[r].abs());
      }

      final double v1 = schottky.k(jInv, m, An) * sqrtOfAbsMu
          - schottky.K(jInv, m, Bn) / sqrtOfAbsMu;
      final double v2 = schottky.k(jInv, m, Bn) / sqrtOfAbsMu
          - schottky.K(jInv, m, An) * sqrtOfAbsMu;

      if (v1 > v2) {
        return v1 / An.dist(Bn);
      }
      else {
        return v2 / An.dist(Bn);
      }
    }

    final double kappa2() {

      double kappa2 = Double.MAX_VALUE;

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          for (int m = 0; m < numGenerators; m++) {
            for (int n = 0; n < numGenerators; n++) {

              if (n != m || i == j) {
                double aKappa = kappaRBar(j, m, i, n);

                if (aKappa < kappa2) {
                  kappa2 = aKappa;
                }
              }
            }
          }
        }
      }
      return kappa2;
    }


    final void H1st(SchottkyGroupElement element) {

      element.diff( B, A, H );

      H.assignDivide( z.re - element.imageOfB[n].re, z.im - element.imageOfB[n].im);
      H.assignDivide( z.re - element.imageOfA[n].re, z.im - element.imageOfA[n].im);

      s.assignPlus( H );
//
//      final double re1 = z.re - element.imageOfB[n].re;
//      final double im1 = z.im - element.imageOfB[n].im;
//
//      final double re2 = z.re - element.imageOfA[n].re;
//      final double im2 = z.im - element.imageOfA[n].im;
//
//      final double as1 = re1 * re1 + im1 * im1;
//      final double as2 = re2 * re2 + im2 * im2;
//
//      s.re += H.re = re1 / as1 - re2 / as2;
//      s.im += H.im = im2 / as2 - im1 / as1;
    }

    final void of1stKind(final SchottkyGroupElement element) {

      if (element.updateID != updateID) {
        schottky.updateElement(element);
      }

      // remark: estimates only valid for word length greater or equal 2,
      // therefore norm is inifinity for generators and id
      final double error = L1[n] * element.norm * rho(element);
 
      if (error * noe[element.wordLength] < acc || error < eps ) {
      	acc += acc / noe[element.wordLength] - error;
      	return;
      }
      
      H1st(element);

      if (element.child == null) {
        schottky.createLeftChilds(element);

      }
      final SchottkyGroupElement[] child = element.child;

      final int numOfChilds = child.length;

      for (int i = 0; i < numOfChilds; i++) {
        of1stKind(child[i]);
      }
    }

    final void of1stKind
        (final Complex r,
         final Complex z,
         final int n, final double accuracy) {

      if (updateID != schottky.updateID)
        update();

      s.assign(0);

      this.z.assign(z);
      this.n = n;
      this.noe = schottky.numOfElementsOfCosetWithWordLength;

      A.assign(schottky.fixpoint[n][0]);
      B.assign(schottky.fixpoint[n][1]);

      H1st(schottky.id);

      if (numGenerators > 1) {

        prepareRho(z);

        this.acc = accuracy;
        this.eps = accuracy / schottky.maxNumOfElements;
        
        for (int i = 0; i < numGenerators; i++) {
          if (i != n) {
            of1stKind(schottky.generator[i]);
            of1stKind(schottky.generatorInv[i]);
          }
        }
      }

      if(this.acc < 0 ) // this test is needed because of the eps crieteria
      	throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
      
      r.assign(s);
    }

    Of1stKindAnalysis of1stKindAnalysis = new Of1stKindAnalysis();

    /**
     * This handles the analysis of the differentials of 1st kind.
     */
    class Of1stKindAnalysis implements Serializable {

      double error, absSeries, waste;
      int maxWordLength, minWordLength, numOfTerms, numOfGiveUps;
      
      boolean adaptAcc = false;

      final void of1stKind(final SchottkyGroupElement element) {

        if (element.updateID != updateID) {
          schottky.updateElement(element);
        }

        // remark: estimates only valid for word length greater or equal 2,
        // therefore norm is inifinity for generators and id
        //final double error = L1[n] * element.norm * rho(element);
        final double error = element.LOfInverse[n] * element.norm * rho(element);

        if ( error * noe[element.wordLength] < acc) {

          waste += acc / noe[element.wordLength] - error;

          if( adaptAcc )
            acc += acc / noe[element.wordLength] - error;

          this.error += error;

          if (element.wordLength < minWordLength)
            minWordLength = element.wordLength;

          if (element.wordLength > maxWordLength)
            maxWordLength = element.wordLength;

          return;
        }

        if ( error < eps ) {
        	acc += acc / noe[element.wordLength] - error;
        	numOfGiveUps++;
        	return;
        }
        
        H1st(element);

        absSeries += H.abs();

        numOfTerms++;

        if (element.child == null) {
          schottky.createLeftChilds(element);

        }
        final SchottkyGroupElement[] child = element.child;

        final int numOfChilds = child.length;

        for (int i = 0; i < numOfChilds; i++) {
          of1stKind(child[i]);
        }
      }

      final void eval
          (final Complex r,
           final Complex z,
           final int n, final double accuracy) {

        eval( r, z, n, accuracy, false );
      }

        final void eval
               (final Complex r,
                final Complex z,
                final int n, final double accuracy, final boolean adaptAcc ) {

        if (updateID != schottky.updateID)
          update();


        AbelianDifferential.this.z.assign(z);
        AbelianDifferential.this.n = n;
        AbelianDifferential.this.noe = schottky.numOfElementsOfCosetWithWordLength;

        this.adaptAcc = adaptAcc; 
        
        s.assign(0);

        A.assign( schottky.fixpoint[n][0] );
        B.assign( schottky.fixpoint[n][1] );

        H1st(schottky.id);

        error = waste = 0;
        absSeries = H.abs();

        numOfTerms = 1;
        numOfGiveUps = 0;
        
        if (numGenerators > 1) {

          prepareRho(z);

          AbelianDifferential.this.z.assign(z);
          AbelianDifferential.this.n = n;
          AbelianDifferential.this.acc = accuracy;
          AbelianDifferential.this.eps = accuracy / schottky.maxNumOfElements;
          
          maxWordLength = 1;
          minWordLength = Integer.MAX_VALUE;

          for (int i = 0; i < numGenerators; i++) {
            if (i != n) {
              of1stKind(schottky.generator[i]);
              of1stKind(schottky.generatorInv[i]);
            }
          }
        } else {
          maxWordLength = 0;
          minWordLength = 0;
        }


        if(AbelianDifferential.this.acc < 0 ) // this test is needed because of the eps crieteria
        	throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
        
        r.assign(AbelianDifferential.this.s);
      }

      final Complex dummy = new Complex();

      final int numOfUsedTerms(final Complex z,
                               final int n, final double accuracy, final boolean adaptAcc) {
        eval(dummy, z, n, accuracy, adaptAcc);
        return numOfTerms;
      }

      final int minWordLength(final Complex z,
                              final int n, final double accuracy, final boolean adaptAcc) {
        eval(dummy, z, n, accuracy, adaptAcc);
        return minWordLength;
      }

      final int maxWordLength(final Complex z,
                              final int n, final double accuracy, final boolean adaptAcc) {
        eval(dummy, z, n, accuracy, adaptAcc);
        return maxWordLength;
      }

      final double error(final Complex z,
                                   final int n, final double accuracy,
                                   final boolean adaptAcc) {
        eval(dummy, z, n, accuracy, adaptAcc);
        return error;
      }

      final double waste(final Complex z,
                                        final int n, final double accuracy, final boolean adaptAcc) {
                       eval(dummy, z, n, accuracy, adaptAcc);
             return waste;
           }

      final double absSeries(final Complex z,
                             final int n, final double accuracy, final boolean adaptAcc) {
        eval(dummy, z, n, accuracy, adaptAcc);
        return absSeries;
      }

    }


    final void H3rd(SchottkyGroupElement element) {

      element.diff(B, A, H);

      element.applyTo(A,a);
      element.applyTo(B,b);

      H.assignDivide(z.re - b.re, z.im - b.im);
      H.assignDivide(z.re - a.re, z.im - a.im);

      s.assignPlus(H);
    }

    final void of3rdKind(final SchottkyGroupElement element) {

      if (element.updateID != updateID) {
        schottky.updateElement(element);
      }

      // remark: estimates only valid for word length greater or equal 2,
      // therefore norm is inifinity for generators and id
      final double error = L1[n] * element.norm * rho(element);

      if (error * noe[element.wordLength] < acc || error < eps ) {
      	acc += acc / noe[element.wordLength] - error;
      	return;
      }
      
      H3rd(element);

      if (element.child == null) {
        schottky.createLeftChilds(element);
      }

      final SchottkyGroupElement[] child = element.child;

      final int numOfChilds = child.length;

      for (int i = 0; i < numOfChilds; i++) {
        of3rdKind(child[i]);
      }
    }

    final void of3rdKind
        (final Complex r,
         final Complex z, final Complex A, final Complex B,
         final double accuracy) {

      if (updateID != schottky.updateID)
        update();

      s.assign(0);

      this.A.assign(A);
      this.B.assign(B);

      this.z.assign(z);
      this.noe = schottky.numOfElementsWithWordLength;

      this.acc = accuracy;
      this.eps = accuracy / schottky.maxNumOfElements;
      
      prepareRho(z);

      of3rdKind( schottky.id );

      if(this.acc < 0 ) // this test is needed because of the eps crieteria
      	throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
      
      r.assign(s);
    }

    final void gamma(final SchottkyGroupElement element) {

      if (element.updateID != updateID) {
        schottky.updateElement(element);
      }

      // remark: estimates only valid for word length greater or equal 2,
      // therefore norm is inifinity for generators and id
      final double error = f2 * element.norm;

      if (error * noe[element.wordLength] < acc || error < eps ) {
      	acc += acc / noe[element.wordLength] - error;
      	return;
      }
      
      element.getInverseOfCSqr(H);

      s.assignPlus(H);

      if (element.child == null) {
        schottky.createLeftChilds(element);

      }
      final SchottkyGroupElement[] child = element.child;

      final int numOfChilds = child.length;

      for (int i = 0; i < numOfChilds; i++) {
        gamma(child[i]);
      }
    }

    final void gamma(final Complex r, final double accuracy) {

      if (updateID != schottky.updateID)
        update();

      s.assign(0);

      this.noe = schottky.numOfElementsWithWordLength;
      this.acc = accuracy;
      this.eps = accuracy / schottky.maxNumOfElements;
      
      for (int i = 0; i < numGenerators; i++) {
        gamma(schottky.generator[i]);
        gamma(schottky.generatorInv[i]);
      }

      if(this.acc < 0 ) // this test is needed because of the eps crieteria
      	throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
      
      r.assign(s);
    }
    
    final void chi(final SchottkyGroupElement element) {
    	
    	if (element.updateID != updateID) {
    		schottky.updateElement(element);
    	}
    	
    	element.getChi( H );
    	
    	if( acc < element.wordLength )
    		return;
    	
    	s.assignPlus(H);
    	
    	if (element.child == null) {
    		schottky.createLeftChilds(element);
    	}
    	
    	final SchottkyGroupElement[] child = element.child;
    	
    	final int numOfChilds = child.length;
    	
    	for (int i = 0; i < numOfChilds; i++) {
    		chi(child[i]);
    	}
    }
    
    final void chi(final Complex r, final double accuracy) {
    	
    	if (updateID != schottky.updateID)
    		update();
    	
    	this.acc = accuracy;
    	this.eps = accuracy / schottky.maxNumOfElements;
    	
    	schottky.id.getChi( s );
    	
    	for (int i = 0; i < numGenerators; i++) {
    		chi(schottky.generator[i]);
    		chi(schottky.generatorInv[i]);
    	}
    	
    	if(this.acc < 0 ) // this test is needed because of the eps crieteria
    		throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
    	
    	r.assign(s);
    }

//
//    final void chi(final SchottkyGroupElement element) {
//
//      if (element.updateID != updateID) {
//        schottky.updateElement(element);
//      }
//
//      // remark: estimates only valid for word length greater or equal 2,
//      // therefore norm is inifinity for generators and id
//      final double error = 2 * rho(element) * element.norm / d2ForZero
//          * ( element.norm / d2ForZero +  maxInIsometricCircles );
//
//      System.out.println( d2ForZero + " " + error );
//
//      if (error * noe[element.wordLength] < acc) {
//        return;
//      }
//
//      schottky.id.getChi( H );
//
//      s.assignPlus(H);
//
//      if (element.child == null) {
//        schottky.createLeftChilds(element);
//      }
//
//      final SchottkyGroupElement[] child = element.child;
//
//      final int numOfChilds = child.length;
//
//      for (int i = 0; i < numOfChilds; i++) {
//        chi(child[i]);
//      }
//    }
//
//    final void chi(final Complex r, final double accuracy) {
//
//      if (updateID != schottky.updateID)
//        update();
//
//      schottky.sigma.update();
//      maxInIsometricCircles = schottky.sigma.maxInIsometricCircles;
//
//      prepareRho( ZERO );
//
//      // compute d2 for zero,
//      if (z.equals( schottky.zOfK2) ) {
//        // default, because we called prepareRho
//        final double[][] k2 = schottky.k2(z);
//        d2ForZero = Double.MIN_VALUE;
//        for (int j = 0; j < 2; j++) {
//          for (int m = 0; m < numGenerators; m++) {
//            if( d2ForZero>k2[j][n] )
//              d2ForZero = k2[j][n];
//          }
//        }
//      } else d2ForZero = schottky.d2( ZERO );
//
//
//      this.noe = schottky.numOfElementsWithWordLength;
//      this.acc = accuracy;
//
//      schottky.id.getChi( s );
//
  //      for (int i = 0; i < numGenerators; i++) {
//        chi(schottky.generator[i]);
//        chi(schottky.generatorInv[i]);
//      }
//
//      r.assign(s);
//    }



}