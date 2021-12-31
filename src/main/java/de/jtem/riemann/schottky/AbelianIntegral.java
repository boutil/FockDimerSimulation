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

class AbelianIntegral {

  final Schottky schottky;

  final double[][] rho;

  final int numGenerators;

  double theta1, q1;

  int updateID;

  AbelianIntegral(Schottky schottky) {
    this.schottky = schottky;
    numGenerators = schottky.numGenerators;
    updateID = schottky.updateID;
    rho = new double[2][numGenerators];
  }

  final Complex zOfRho
      = new Complex(Double.NaN);

  final Complex z = new Complex();
  final Complex p = new Complex();
  final Complex H = new Complex();

  final Complex A = new Complex();
  final Complex B = new Complex();

  final Complex a = new Complex();
  final Complex b = new Complex();

  final Complex d = new Complex();

  int n;

  long [] noe;

  double acc;
  double eps; // = acc / maxNumberOfElements;
  
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

  final void of1stKind(final SchottkyGroupElement element) {

    if (element.updateID != updateID) {
      schottky.updateElement(element);
    }

    element.diff(B, A, d);

    final double error = (Math.abs(d.re) + Math.abs(d.im)) * rho(element);

    if (error * noe[element.wordLength] < acc || error < eps ) {
      acc += acc / noe[element.wordLength] - error;
      return;
    }
    
    H.assignDivide(z.re - element.imageOfB[n].re,
                   z.im - element.imageOfB[n].im,
                   z.re - element.imageOfA[n].re,
                   z.im - element.imageOfA[n].im);

    p.assignTimes(H);

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

    p.assignDivide(z.re - schottky.fixpoint[n][1].re,
                   z.im - schottky.fixpoint[n][1].im,
                   z.re - schottky.fixpoint[n][0].re,
                   z.im - schottky.fixpoint[n][0].im);

    // p.assignTimes( -1 ); // exp( pi i )

    if (numGenerators > 1) {

      prepareRho(z);

      this.noe = schottky.numOfElementsOfCosetWithWordLength;

      this.A.assign(schottky.fixpoint[n][0]);
      this.B.assign(schottky.fixpoint[n][1]);

      this.z.assign(z);
      this.n = n;
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
    
    r.assignLog(p);
  }

  Of1stKindAnalysis of1stKindAnalysis = new Of1stKindAnalysis();

  /**
   * This handles the analysis of the integrals of 1st kind.
   */
  class Of1stKindAnalysis
      implements Serializable {

    double error, absSeries;
    int maxWordLength, minWordLength, numOfTerms, numOfGiveUps;
    
    Complex logOfH = new Complex();

    final void of1stKind(final SchottkyGroupElement element) {

      if (element.updateID != updateID) {
        schottky.updateElement(element);
      }

      element.diff(B, A, d);

      final double error = ( Math.abs(d.re) + Math.abs(d.im) ) * rho(element);

      if ( error * noe[element.wordLength] < acc) {

        acc += acc / noe[element.wordLength] - error;

        this.error += error;

        if (element.wordLength < minWordLength) {
          minWordLength = element.wordLength;
        }
        if (element.wordLength < maxWordLength) {
          maxWordLength = element.wordLength;
        }

        return;
      }

      if ( error < eps ) {
      	acc += acc / noe[element.wordLength] - error;
      	numOfGiveUps++;
      	return;
      }
      
      H.assignDivide(z.re - element.imageOfB[n].re,
                     z.im - element.imageOfB[n].im,
                     z.re - element.imageOfA[n].re,
                     z.im - element.imageOfA[n].im);

      p.assignTimes(H);

      logOfH.assignLog( H );

      absSeries += logOfH.abs();

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

      p.assignDivide(z.re - schottky.fixpoint[n][1].re,
                     z.im - schottky.fixpoint[n][1].im,
                     z.re - schottky.fixpoint[n][0].re,
                     z.im - schottky.fixpoint[n][0].im);

      //p.assignTimes( -1 ); // exp( pi i )

      error = 0;

      numOfTerms = 1;
      numOfGiveUps =0;
      
      logOfH.assignLog( p );

      absSeries = logOfH.abs(); // + Math.PI;


      if (numGenerators > 1) {

        minWordLength = Integer.MAX_VALUE;

        prepareRho(z);

        AbelianIntegral.this.noe = schottky.numOfElementsOfCosetWithWordLength;

        AbelianIntegral.this.A.assign(schottky.fixpoint[n][0]);
        AbelianIntegral.this.B.assign(schottky.fixpoint[n][1]);

        AbelianIntegral.this.z.assign(z);
        AbelianIntegral.this.n = n;
        AbelianIntegral.this.acc = accuracy;
        AbelianIntegral.this.eps = accuracy / schottky.maxNumOfElements;
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

      if(AbelianIntegral.this.acc < 0 ) // this test is needed because of the eps crieteria
      	throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
      
      r.assignLog(p);
    }

    final Complex dummy = new Complex();

    final int numOfUsedTerms(final Complex z,
                             final int n, final double accuracy) {
      eval(dummy, z, n, accuracy);
      return numOfTerms;
    }

    final int minWordLength(final Complex z,
                            final int n, final double accuracy) {
      eval(dummy, z, n, accuracy);
      return minWordLength;
    }

    final int maxWordLength(final Complex z,
                            final int n, final double accuracy) {
      eval(dummy, z, n, accuracy);
      return maxWordLength;
    }

    final double error(final Complex z,
                       final int n, final double accuracy) {
      eval(dummy, z, n, accuracy);
      return error;
    }

    final double absSeries(final Complex z,
                           final int n, final double accuracy) {
      eval(dummy, z, n, accuracy);
      return absSeries;
    }

  }


  final void of3rdKind(final SchottkyGroupElement element) {

    if (element.updateID != updateID) {
      schottky.updateElement(element);
    }

    element.diff(B, A, d);

    final double error = (Math.abs(d.re) + Math.abs(d.im)) * rho(element);

    if (error * noe[element.wordLength] < acc || error < eps ) {
    	acc += acc / noe[element.wordLength] - error;
    	return;
    }
    
    element.applyTo(A, a);
    element.applyTo(B, b);

    H.assignDivide(z.re - b.re, z.im - b.im,
                   z.re - a.re, z.im - a.im);

    p.assignTimes(H);

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
       final Complex z,
       final Complex A, final Complex B,
       final double accuracy) {

    this.noe = schottky.numOfElementsWithWordLength;

    this.A.assign(A);
    this.B.assign(B);

    this.z.assign(z);
    this.acc = accuracy;
    this.eps = accuracy / schottky.maxNumOfElements;
    
    prepareRho(z);

    p.assign(1);

    of3rdKind(schottky.id);
    
    if(this.acc < 0 ) // this test is needed because of the eps crieteria
    	throw new RuntimeException( "could not evaluate series because of numerical instabilities" );
    
    r.assignLog(p);
  }


}