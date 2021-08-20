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

class Sigma
    implements Serializable {

  final Schottky schottky;

  final int numGenerators;

  double theta1, q1, maxInIsometricCircles;

  int updateID;

  Sigma(Schottky schottky) {
    this.schottky = schottky;
    numGenerators = schottky.numGenerators;
  }

  void update() {

    if (updateID == schottky.updateID) {
      return;
    }

    updateID = schottky.updateID;

    theta1 = schottky.theta1;
    q1 = schottky.q1;

    maxInIsometricCircles = maxInIsometricCircles();
  }

  /**
   * Computes biggest value in isometric circles.
   * @return biggest value in isometric circles.
   */
  double maxInIsometricCircles() {
    double max=0;
    for( int i=0; i<2; i++ )
      for( int n=0; n<numGenerators; n++ ) {
        final double value = schottky.radius[n] + schottky.center[n][i].abs();
        if( max < value )
          max = value;
      }
    return max;
  }

  int n; // index of coset
  int k; // power

  long[] noe;

  final Complex s = new Complex();
  final Complex d = new Complex();
  final Complex z = new Complex();
  final Complex w = new Complex();

  double acc;
  double factor;

  final private void sigma(SchottkyGroupElement element) {

    if (element.updateID != updateID) {
      schottky.updateElement(element);
    }

    element.diff(z, w, d);

    
    final double error = factor * noe[element.wordLength]
        * (Math.abs(d.re) + Math.abs(d.im));

    if( element.wordLength > 50 ) {
    	System.out.println( "factor=" + factor);
    	System.out.println( "z=" + z);
    	System.out.println( "d=" + d);
    	System.out.println("error="+error);
    	
    	throw new RuntimeException("ups");
    }
    if (error < acc) {
      return;
    }

    s.assignPlus(d);

    if (element.child == null) {
      schottky.createLeftChilds(element);

    }
    final SchottkyGroupElement[] child = element.child;

    final int numOfChilds = child.length;

    for (int i = 0; i < numOfChilds; i++) {
      sigma(child[i]);
    }
  }

  final void eval(Complex r, int n, double accuracy) {

    this.n = n;
    this.noe = schottky.numOfElementsOfCosetWithWordLength;
    this.acc = accuracy;
    this.factor = 1 / (1 - q1);

    this.z.assign( schottky.fixpoint[n][0] );
    this.w.assign( schottky.fixpoint[n][1] );

    s.assignMinus(z, w);

    for (int i = 0; i < numGenerators; i++) {
      if (i != n) {
        sigma(schottky.generator[i]);
        sigma(schottky.generatorInv[i]);
      }
    }

    r.assign(s);
  }

  final void eval( Complex r,
                   Complex z, Complex w, double accuracy) {

    this.noe = schottky.numOfElementsWithWordLength;
    this.acc = accuracy;
    this.factor = 1 / (1 - q1);

    this.z.assign(z);
    this.w.assign(w);

    s.assignMinus(z, w);

    for (int i = 0; i < numGenerators; i++) {
      sigma(schottky.generator[i]);
      sigma(schottky.generatorInv[i]);
    }

    r.assign(s);
  }

  Analysis analysis = new Analysis();

  /**
   * This handles the analysis of the differentials of 1st kind.
   */
  class Analysis
      implements Serializable {

    double error, absSeries;
    int maxWordLength, minWordLength, numOfTerms, numOfSmallTerms;

    final void sigma(final SchottkyGroupElement element) {

      if (element.updateID != updateID) {
        schottky.updateElement(element);
      }

      element.diff(z, w, d);

      final double error = factor * noe[element.wordLength]
          * (Math.abs(d.re) + Math.abs(d.im));

      if (error < acc) {

        this.error += factor * (Math.abs(d.re) + Math.abs(d.im));

        if (element.wordLength < minWordLength) {
          minWordLength = element.wordLength;

        }
        if (element.wordLength > maxWordLength) {
          maxWordLength = element.wordLength;

        }
        return;
      }

      if (d.re + s.re == s.re && d.im + s.im == s.im) {
        numOfSmallTerms++;

      }
      s.assignPlus(d);

      absSeries += d.abs();

      numOfTerms++;

      if (element.child == null) {
        schottky.createLeftChilds(element);

      }
      final SchottkyGroupElement[] child = element.child;

      final int numOfChilds = child.length;

      for (int i = 0; i < numOfChilds; i++) {
        sigma(child[i]);
      }
    }

    final void eval(Complex r,
                    Complex z, Complex w, double accuracy) {

      Sigma.this.n = n;
      Sigma.this.noe = schottky.numOfElementsWithWordLength;
      Sigma.this.acc = accuracy;
      Sigma.this.factor = 1 / (1 - q1);

      Sigma.this.z.assign(z);
      Sigma.this.w.assign(w);

      s.assignMinus(z, w);

      error = 0;
      absSeries = s.abs();

      numOfTerms = 1;
      numOfSmallTerms = 0;
      maxWordLength = 1;
      minWordLength = Integer.MAX_VALUE;

      for (int i = 0; i < numGenerators; i++) {
        sigma(schottky.generator[i]);
        sigma(schottky.generatorInv[i]);
      }

      r.assign(s);
    }

    final Complex dummy = new Complex();

    final int numOfUsedTerms(final Complex z, final Complex w,
                             final double accuracy) {
      eval(dummy, z, w, accuracy);
      return numOfTerms;
    }

    final int numOfSmallTerms(final Complex z, final Complex w,
                              final double accuracy) {
      eval(dummy, z, w, accuracy);
      return numOfSmallTerms;
    }

    final int minWordLength(final Complex z, final Complex w,
                            final double accuracy) {
      eval(dummy, z, w, accuracy);
      return minWordLength;
    }

    final int maxWordLength(final Complex z, final Complex w,
                            final double accuracy) {
      eval(dummy, z, w, accuracy);
      return maxWordLength;
    }

    final double error(final Complex z, final Complex w,
                       final double accuracy) {
      eval(dummy, z, w, accuracy);
      return error;
    }

    final double absSeries(final Complex z, final Complex w,
                           final double accuracy) {
      eval(dummy, z, w, accuracy);
      return absSeries;
    }

  }

  final private void sigmaPow(SchottkyGroupElement element) {

    if (element.updateID != updateID) {
      schottky.updateElement(element);
    }

    element.diff(z, w, d);

    final double error = factor * noe[element.wordLength]
        * (Math.abs(d.re) + Math.abs(d.im));

    if (error < acc) {
      return;
    }

    element.diffPow( z, w, k, d, d );

    s.assignPlus(d);

    if (element.child == null) {
      schottky.createLeftChilds(element);

    }
    final SchottkyGroupElement[] child = element.child;

    final int numOfChilds = child.length;

    for (int i = 0; i < numOfChilds; i++) {
      sigmaPow(child[i]);
    }
  }

  final void evalPow(Complex r, int n, int k, double accuracy) {

    this.n = n;
    this.noe = schottky.numOfElementsOfCosetWithWordLength;
    this.acc = accuracy;
    this.factor = k * maxInIsometricCircles / (1 - q1);

    this.z.assign( schottky.fixpoint[n][0]);
    this.w.assign( schottky.fixpoint[n][1]);
    this.k = k;

    schottky.id.diffPow(z, w, k, s);

    for (int i = 0; i < numGenerators; i++) {
      if (i != n) {
        sigmaPow(schottky.generator[i]);
        sigmaPow(schottky.generatorInv[i]);
      }
    }

    r.assign(s);
  }

  final void eval(Complex r,
                  Complex z, Complex w, int k, double accuracy) {

    this.noe = schottky.numOfElementsWithWordLength;
    this.acc = accuracy;
    this.factor = k * maxInIsometricCircles / (1 - q1);

    this.z.assign(z);
    this.w.assign(w);
    this.k = k;

    schottky.id.diffPow(z, w, k, s);

    for (int i = 0; i < numGenerators; i++) {
      sigmaPow(schottky.generator[i]);
      sigmaPow(schottky.generatorInv[i]);
    }

    r.assign(s);
  }

}