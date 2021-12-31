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

import de.jtem.mfc.field.Complex;

/**
 * @version 0.5
 * @author schmies
 *
 */
public class SchottkyAnalysis
    extends Schottky {

  private static final long serialVersionUID = 1L;

  /**
   *
   */
  public SchottkyAnalysis() {
    this(2);
  }

  /**
   * @param genus
   */
  public SchottkyAnalysis(int genus) {
    this(getDefaultUniformizationData(genus));
  }

  /**
   * @param uniformizationData
   */
  public SchottkyAnalysis(double[] uniformizationData) {
    this(uniformizationData, 1e-7);
  }

  /**
   * @param uniformizationData
   * @param accuracy
   */
  public SchottkyAnalysis(double[] uniformizationData, double accuracy) {
    super(uniformizationData, accuracy);
  }

  public boolean getUseFancyError() {
    return useFancyError;
  }

  public void setUseFancyError(boolean v) {
    useFancyError = v;
  }

  public final double getTheta1() {
    return theta1;
  }

  public final double getQ1() {
    return q1;
  }

  public final double getKappa2() {
    return abelianDifferential.getKappa2();
  }

  public final double getQ2() {
    return abelianDifferential.getQ2();
  }

  /**
   * Returns the minimal distance of the images of the fundamental
   * domain F under all transformation of word length 2.
   */
  final double getD2() {

    double d2 = Double.MAX_VALUE;

    for (int n = 0; n < numGenerators; n++) {
      for (int i = 0; i < 2; i++) {

        Complex[] someInnerCenter = innerCenter[n][i];
        double[] someInnerRadius = innerRadius[n][i];

        for (int k = 0; k < 2 * numGenerators - 1; k++) {

          double tmpD = radius[n]
              - someInnerCenter[k].dist(center[n][i])
              - someInnerRadius[k];

          if (tmpD < d2) {
            d2 = tmpD;
          }
        }
      }
    }
    return d2;
  }

  public double[] getL1() {
    return abelianDifferential.L1;
  }

  public double d(SchottkyGroupElement sigma) {

    double dist = Double.MAX_VALUE;

    int n = sigma.left;
    int i = sigma.leftIsInvert == 1 ? 0 : 1;

    for (int m = 0; m < numGenerators; m++) {
      for (int j = 0; j < 2; j++) {

        if (m != n || i == j) {

          double aRadius = sigma.getRadiusOfMappedCircle(center[m][j], radius[m],
              tmp);

          double aDist = radius[n] - center[n][i].dist(tmp) - aRadius;

          // check that it is in the target circle
          if (aDist < 0) {
            System.out.println("computeInnerCircles: fail target");
          }

          if (aDist < dist) {
            dist = aDist;
          }
        }
      }
    }

    return dist;
  }



    final private double eval_d(SchottkyGroupElement element,
                                double d, int wordLength) {

      if (element.wordLength > wordLength) {
        return d;
      }

      if (element.updateID != updateID) {
        updateElement(element);
      }

      if (element.wordLength == wordLength) {
        d = Math.min( d(element), d );
      }

      if (element.child == null) {
        createLeftChilds(element);
      }

      final SchottkyGroupElement[] child = element.child;

      final int numOfChilds = child.length;

      for (int i = 0; i < numOfChilds; i++) {
        d = eval_d(child[i], d, wordLength);
      }
      return d;
    }

    /**
     * Computes d for given word lengths.
     */
    public final double d( int wordLength) {

      if (wordLength < 2) {
        throw new IllegalArgumentException
            ("kappa only defined for word legnth greater then one");
      }

      if (updateID != updateID) {
        update();

      }
      return eval_d( id, Double.MAX_VALUE, wordLength);
    }

    /**
     * Returns group element with provided word.
     * @param word letters a,b,c, ... prescripe generators A,B,C, ... their inverses.
     * @return group element with provided word
     */
    public final SchottkyGroupElement getGroupElement( String word ) {

      char [] letter = word.toCharArray();

      SchottkyGroupElement element = id;

      for( int i=letter.length - 1; i>=0; i-- ) {

        if (element.updateID != updateID) {
          updateElement(element);
        }

        if (element.child == null) {
          createLeftChilds(element);
        }

        element = element.childWithFirstLetter(letter[i]);
      }

      if (element.updateID != updateID) {
        updateElement(element);
      }

      return element;
    }


  public static String inverse( String word ) {

    char [] letter = word.toCharArray();

    char [] inverse = new char[letter.length];

    for (int i = letter.length - 1, j=0; j < letter.length; j++, i--) {

      if( letter[i] <'a' )
        inverse[j] = (char)((int)'a' + letter[i] - 'A');
      else
        inverse[j] = (char)((int)'A' + letter[i] - 'a');
    }

    return new String( inverse );
  }
  }