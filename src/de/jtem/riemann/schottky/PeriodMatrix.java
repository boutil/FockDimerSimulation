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

import de.jtem.blas.ComplexMatrix;
import de.jtem.mfc.field.Complex;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * @author Markus Schmies
 * @version 1.0
 */

  class PeriodMatrix implements Serializable {

    final Schottky schottky;

    final int numGenerators;

    double theta1, q1;

    int updateID;

    PeriodMatrix(Schottky schottky) {
      this.schottky = schottky;
      numGenerators = schottky.numGenerators;
    }

    int n;
    int m;

    double acc;

    double factorForLeftIsNotM;
    double factorForLeftIsM;

    final Complex b = new Complex();
    final Complex p = new Complex();
    final Complex d = new Complex();

    final Complex An = new Complex();
    final Complex Bn = new Complex();
    final Complex Am = new Complex();
    final Complex Bm = new Complex();

    void update() {

      if (updateID == schottky.updateID) {
        return;
      }

      updateID = schottky.updateID;

      theta1 = schottky.theta1;
      q1 = schottky.q1;
    }

    final private void B(final SchottkyGroupElement element) {

      if (element.updateID != updateID) {
        schottky.updateElement(element);
      }

      final long noe = schottky.getNumOfElementsOfCosetWithWordLength(element.wordLength);

      element.diff(Bn, An, d);

      final double dist = Math.abs(d.re) + Math.abs(d.im);
//
//      if( element.left != m )
//        System.out.println( dist + " " + (factorForLeftIsNotM + " " + noe) );
//      else
//        System.out.println( dist + " " + (factorForLeftIsM + "  " + noe) );

      if (element.left != m && factorForLeftIsNotM * dist * noe < acc ||
          element.left == m && factorForLeftIsM * dist * noe < acc) {
        return;
      }

      if (element.left != m) {
        p.assignCrossRatio( Am, element.imageOfB[n],
                            Bm, element.imageOfA[n]);

        b.assignTimes(p);
      }

      if (element.child == null) {
        schottky.createLeftChilds(element);

      }
      final SchottkyGroupElement[] child = element.child;

      final int numOfChilds = child.length;

      for (int i = 0; i < numOfChilds; i++) {
        B(child[i]);
      }
    }

    final double k1(final int m, final Complex P) {

      double min = Double.MAX_VALUE;

      for (int n = 0; n < numGenerators; n++) {
        if (n != m) {
          final double k1n = Math.min(schottky.k(0, n, P), schottky.k(1, n, P));

          if (min > k1n) {
            min = k1n;
          }
        }
      }

      return min;
    }

    /**
     * Computes period matrix with prescribed accuracy.
     * @param B period matrix on output
     * @param accuracy of computed period matrix
     */
    void eval(ComplexMatrix B, double accuracy) {

      B.newSize(numGenerators);

      if (numGenerators == 1) { // treat genus one case seperatly
        b.assignLog(schottky.mu[0]);
        B.set(0, 0, b);
        return;
      }

      this.acc = accuracy;

      final double v = theta1 * theta1;

      final double rOfV = schottky.r(v);
      final double sumForLeftIsM = (2 * numGenerators - 2) * rOfV;
      final double sumForLeftIsNotM = (2 * numGenerators - 4) * rOfV +
          schottky.rPlus(rOfV, v) +
          schottky.rMinus(rOfV, v);

      for (m = 0; m < numGenerators; m++) {

        Am.assign(schottky.fixpoint[m][0]);
        Bm.assign(schottky.fixpoint[m][1]);

        factorForLeftIsM = factorForLeftIsNotM
            = 1 / k1(m, Am ) + 1 / k1(m, Bm);

        factorForLeftIsM *= sumForLeftIsM;
        factorForLeftIsNotM *= sumForLeftIsNotM;

        for (n = 0; n < numGenerators; n++) {

          if (m > n) {
            B.get(m, n, b);
          }
          else {
            An.assign(schottky.fixpoint[n][0]);
            Bn.assign(schottky.fixpoint[n][1]);

            if (m < n) {
              b.assignCrossRatio( Am, Bn, Bm, An );
            }
            else { // n==m
              b.assign(schottky.mu[n]);
            }

            for (int i = 0; i < numGenerators; i++) {
              if (i != n) {
                B(schottky.generator[i]);
                B(schottky.generatorInv[i]);
              }
            }

            b.assignLog();
          }

          B.set(n, m, b);
        }
      }
    }
  }
