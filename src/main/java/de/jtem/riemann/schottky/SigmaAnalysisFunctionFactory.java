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

public class SigmaAnalysisFunctionFactory {

	static final Complex ZERO = new Complex();


  /**
   * Creates function which evaluates numerical the absolute error
   * we make when we use the approxiamation with a prescribed
   * error given by testAcc. This is simply done by measering the distance
   * of the approximation, with prescribed testAcc, to a second
   * approximation, with significant higher accuracy, the controlAcc
   * which is used for the approximation of sigma.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function counting the number of elment in section.
   */
  public static SchottkyAnalysisFunction
      createNumericalErrorFunction
      (Schottky schottky, double testAcc, double controlAcc ) {

    return new NumericalErrorFunction(schottky, testAcc, controlAcc, ZERO);
  }
  
  private static final class NumericalErrorFunction extends SchottkyAnalysisFunction {
		  private final Complex ZERO;
		  private NumericalErrorFunction(Schottky schottky, double testAcc, double controlAcc, Complex ZERO) {
			  super(schottky, testAcc, controlAcc);
			  this.ZERO = ZERO;
		  }
		  public double getControlAcc() {
			  return super.getControlAcc();
			}
		  public void setControlAcc(double controlAcc) {
			  super.setControlAcc(controlAcc);
			}
		  public double getTestAcc() {
			  return super.getTestAcc();
			}
		  public void setTestAcc(double testAcc) {
			  super.setTestAcc(testAcc);
			}
		  public double eval(Complex z) {
		
			  if (controlAcc < schottky.getAccuracy()) {
				throw new IllegalArgumentException
					("controlAcc exeeds accuracy of analysing schottky instance");
			  }
		
			  schottky.sigma( controlValue, z, ZERO, controlAcc);
			  schottky.sigma( testValue,    z, ZERO, acc);
		
			  return controlValue.dist(testValue);
			}
	  }


  /**
   * Creates function counting the number of elements in subset which
   * is used for the approximation of sigma.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function counting the number of elment in used subset.
   */
  public static SchottkyAnalysisFunction
      createNumOfUsedTermsFunction
      (Schottky schottky, double testAcc) {

    return new NumOfUsedTermsFunction(schottky, testAcc, ZERO);
  }

  private static final class NumOfUsedTermsFunction extends SchottkyAnalysisFunction.ForProperty {
	  private final Complex ZERO;
	  private NumOfUsedTermsFunction(Schottky schottky, double testAcc, Complex ZERO) {
		  super(schottky, testAcc);
		  this.ZERO = ZERO;
	  }
	  public double eval(Complex z) {
		  return schottky.sigma.analysis.numOfUsedTerms(z, ZERO, acc);
		}
  }

  /**
   * Creates function counting the number of elements in subset which
   * deliver a term below the numerical resolution.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function counting the number of small elments in used subset.
   */
  public static SchottkyAnalysisFunction
      createNumOfSmallTermsFunction
      (Schottky schottky, double testAcc) {

    return new NumOfSmallTermsFunction(schottky, testAcc, ZERO);
      }


	private static final class NumOfSmallTermsFunction extends SchottkyAnalysisFunction.ForProperty {
		private final Complex ZERO;
		private NumOfSmallTermsFunction(Schottky schottky, double testAcc, Complex ZERO) {
			super(schottky, testAcc);
			this.ZERO = ZERO;
		}
		public double eval(Complex z) {
			return schottky.sigma.analysis.numOfSmallTerms(z, ZERO, acc);
		  }
	}


  /**
   * Creates function which delivers the smallest word length of all
   * element in section of the the subset which
   * is used for the approximation sigma.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function delivering smallest word length of section
   */
  public static SchottkyAnalysisFunction
      createMinWordLengthFunction
      (Schottky schottky, double testAcc) {

    return new MinWordLengthFunction(schottky, testAcc, ZERO);
  }


  private static final class MinWordLengthFunction extends SchottkyAnalysisFunction.ForProperty {
	  private final Complex ZERO;
	  private MinWordLengthFunction(Schottky schottky, double testAcc, Complex ZERO) {
		  super(schottky, testAcc);
		  this.ZERO = ZERO;
	  }
	  public double eval(Complex z) {
		  return schottky.sigma.analysis.minWordLength(z,ZERO,acc);
		}
  }

  /**
   * Creates function which checks the analysis object for
   * sigma.
   * This is done by comparing the result of the regulare
   * routine with the one for the analysis. This values must
   * coinside, otherwise something is wrong.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function delivering smallest word length of section
   */
  public static SchottkyAnalysisFunction
      createAnalysisCheckFunction
      (Schottky schottky, double testAcc) {

    return new AnalysisCheckFunction(schottky, testAcc, ZERO);
  }

  private static final class AnalysisCheckFunction extends de.jtem.riemann.schottky.SchottkyAnalysisFunction.ForProperty {
	  private final Complex ZERO;
	  final Complex r = new Complex();
	  final Complex a = new Complex();
	  private AnalysisCheckFunction(Schottky schottky, double testAcc, Complex ZERO) {
		  super(schottky, testAcc);
		  this.ZERO = ZERO;
	  }
	  public double eval(Complex z) {
		  schottky.sigma.analysis.eval(a, z, ZERO, acc);
		  schottky.sigma(r, z, ZERO, acc);
		  return a.dist(r);
		}
  }

  /**
   * Creates function which delivers the largest word length of all
   * element in section of the the subset which
   * is used for the approximation of sigma.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function delivering largestest word length of section
   */
  public static SchottkyAnalysisFunction
      createMaxWordLengthFunction
      (Schottky schottky, double testAcc) {

    return new MaxWordLengthFunction(schottky, testAcc, ZERO);
  }

  private static final class MaxWordLengthFunction extends SchottkyAnalysisFunction.ForProperty {
	  private final Complex ZERO;
	  private MaxWordLengthFunction(Schottky schottky, double testAcc, Complex ZERO) {
		  super(schottky, testAcc);
		  this.ZERO = ZERO;
	  }
	  public double eval(Complex z) {
		  return schottky.sigma.analysis.
			  maxWordLength(z, ZERO, acc);
		}
  }

  /**
   * Creates function which delivers the absolute series
   * related to sigma.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function delivering the absolute series
   */
  public static SchottkyAnalysisFunction
      createAbsSeriesFunction
      (Schottky schottky, double testAcc) {

    return new AbsSeriesFunction(schottky, testAcc, ZERO);
  }


  private static final class AbsSeriesFunction extends SchottkyAnalysisFunction.ForProperty {
	  private final Complex ZERO;
	  private AbsSeriesFunction(Schottky schottky, double testAcc, Complex ZERO) {
		  super(schottky, testAcc);
		  this.ZERO = ZERO;
	  }
	  public double eval(Complex z) {
		  return schottky.sigma.analysis.
			  absSeries(z, ZERO, acc);
		}
  }


  /**
   * Creates function which delivers the theoretical error
   * of the approximation of sigma.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function counting the number of elment in section.
   */
  public static SchottkyAnalysisFunction
      createErrorFunction
      (Schottky schottky, double testAcc) {

    return new ErrorFunction(schottky, testAcc, ZERO);
  }

  private static final class ErrorFunction extends SchottkyAnalysisFunction.ForProperty {
	  private final Complex ZERO;
	  private ErrorFunction(Schottky schottky, double testAcc, Complex ZERO) {
		  super(schottky, testAcc);
		  this.ZERO = ZERO;
	  }
	  public double eval(Complex z) {
		  return schottky.sigma.analysis.
			  error(z, ZERO, acc);
		}
  }


    /**
     * Creates function evaluating the absolute value of sigma.
     * @param schottky instance of analysis
     * @param testAcc accuracy of approximation
     * @return function computing the absolute value of sigma.
     */
    public static SchottkyAnalysisFunction
        createAbsValueFunction
        (Schottky schottky, double testAcc ) {

      return new AbsValueFunction(schottky, testAcc, ZERO);
  }

  private static final class AbsValueFunction extends SchottkyAnalysisFunction.ForProperty {
	  private final Complex ZERO;
	  private AbsValueFunction(Schottky schottky, double testAcc, Complex ZERO) {
		  super(schottky, testAcc);
		  this.ZERO = ZERO;
	  }
	  public double eval(Complex z) {
		schottky.sigma(testValue, z, ZERO, acc);
		return testValue.abs();
	  }
  }


  /**
   * Creates function evaluating the absolute value of sigma.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @return function computing the absolute value of sigma.
   */
  public static SchottkyAnalysisFunction
      createAbsValueFunctionForPower
      (Schottky schottky, int pow, double testAcc ) {

    return new AbsValueFunctionForPower(schottky, testAcc, pow, ZERO);
  }

  private static final class AbsValueFunctionForPower extends SchottkyAnalysisFunction.ForIndexedProperty {
		private final Complex ZERO;
		private AbsValueFunctionForPower(Schottky schottky, double testAcc, int n, Complex ZERO) {
			super(schottky, testAcc, n);
			this.ZERO = ZERO;
		}
		public double eval(Complex z) {
			schottky.sigma.eval(testValue, z, ZERO, n, acc);
			return testValue.abs();
		  }
	}

}
