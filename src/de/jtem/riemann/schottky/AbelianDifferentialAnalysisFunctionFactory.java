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

public class AbelianDifferentialAnalysisFunctionFactory {

  private final static class ArgValueFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		private ArgValueFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public double eval(Complex z) {
          schottky.abelianDifferentialOf1stKind(testValue, z, n, acc);
          return testValue.arg();
        }
	}

	private final static class BenchMarkFunction0 extends SchottkyAnalysisFunction.BenchMarkFunction {
		private BenchMarkFunction0(SchottkyAnalysisFunction function) {
			super(function);
		}
		public double eval(Complex z) {
          if (function.acc < schottky.getAccuracy()) {
            throw new IllegalArgumentException
                ("testAcc exeeds accuracy of analysing schottky instance");
          }

          long maxDuration = 0, totalDuration = 0;

          for (int i = 0; i < 11; i++) {

            long timeMillis = System.currentTimeMillis();

            for (int j = 0; j < 10; j++) {
              function.eval(z);
            }

            long duration = System.currentTimeMillis() - timeMillis;

            totalDuration += duration;

            if (duration > maxDuration) {
              maxDuration = duration;
            }
          }

          return (totalDuration - maxDuration) / 100.;
        }
	}

	private final static class AbsValueFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		private AbsValueFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public double eval(Complex z) {
		    schottky.abelianDifferentialOf1stKind(testValue, z, n, acc);
		    return testValue.abs();
		  }
	}

	private final static class WasteFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		private WasteFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public double eval(Complex z) {
            return schottky.abelianDifferential.of1stKindAnalysis.
                waste(z, n, acc, false);
          }
	}

	private final static class ErrorFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		boolean adaptAcc = false;
		private ErrorFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public boolean isAdaptAcc() {
		  return adaptAcc;
		}
		public void setAdaptAcc( boolean adaptAcc ) {
		  this.adaptAcc = adaptAcc;
		}
		public double eval(Complex z) {
          return schottky.abelianDifferential.of1stKindAnalysis.
              error(z, n, acc, adaptAcc);
        }
	}

	private final static class AbsSeriesFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		private AbsSeriesFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public double eval(Complex z) {
		    return schottky.abelianDifferential.of1stKindAnalysis.
		        absSeries(z, n, acc, false);
		  }
	}

	private final static class MaxWordLengthFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		private MaxWordLengthFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public double eval(Complex z) {
		    return schottky.abelianDifferential.of1stKindAnalysis.
		        maxWordLength(z, n, acc, false);
		  }
	}

	private final static class AnalysisCheckFunctionFor1stKind extends de.jtem.riemann.schottky.SchottkyAnalysisFunction.ForIndexedProperty {
		final Complex r = new Complex();
		final Complex a = new Complex();
		private AnalysisCheckFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public double eval(Complex z) {
		    schottky.abelianDifferential.of1stKindAnalysis.
		        eval(a, z, n, acc);
		    schottky.abelianDifferentialOf1stKind(r,z,n,acc);
		    return a.dist(r);
		  }
	}

	private final static class MinWordLengthFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		private MinWordLengthFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public double eval(Complex z) {
		    return schottky.abelianDifferential.of1stKindAnalysis.
		        maxWordLength(z, n, acc, false);
		  }
	}

	private final static class NumOfUsedTermsFunctionFor1stKind extends SchottkyAnalysisFunction.ForIndexedProperty {
		boolean adaptAcc = false;
		private NumOfUsedTermsFunctionFor1stKind(Schottky schottky, double testAcc, int n) {
			super(schottky, testAcc, n);
		}
		public boolean isAdaptAcc() {
		    return adaptAcc;
		  }
		public void setAdaptAcc(boolean adaptAcc) {
		    this.adaptAcc = adaptAcc;
		          }
		public double eval(Complex z) {
		    return schottky.abelianDifferential.of1stKindAnalysis.
		        numOfUsedTerms(z, n, acc, adaptAcc);
		  }
	}

	private final static class NumericalErrorOfAbsFunctionFor1stKind extends SchottkyAnalysisFunction {
		boolean adaptAcc = false;
		private NumericalErrorOfAbsFunctionFor1stKind(Schottky schottky, double testAcc, double controlAcc, int n) {
			super(schottky, testAcc, controlAcc, n);
		}
		public boolean isAdaptAcc() {
		     return adaptAcc;
		   }
		public void setAdaptAcc(boolean adaptAcc) {
		     this.adaptAcc = adaptAcc;
		   }
		public double getControlAcc() {
		     return super.getControlAcc();
		   }
		public void setControlAcc(double controlAcc) {
		     super.setControlAcc(controlAcc);
		   }
		public int getN() {
		     return super.getN();
		   }
		public void setN(int n) {
		     super.setN(n);
		   }
		public double getTestAcc() {
		     return super.getTestAcc();
		   }
		public void setTestAcc(double testAcc) {
		     super.setTestAcc(testAcc);
		   }
		public double eval(Complex z) {
		
		     return schottky.abelianDifferential.of1stKindAnalysis.absSeries( z, n, controlAcc, adaptAcc)-
		         schottky.abelianDifferential.of1stKindAnalysis.absSeries( z, n, acc, adaptAcc );
		    }
	}

	private final static class NumErr1stKind extends SchottkyAnalysisFunction {
		private NumErr1stKind(Schottky schottky, double testAcc, double controlAcc, int n) {
			super(schottky, testAcc, controlAcc, n);
		}
		public double getControlAcc() {
		    return super.getControlAcc();
		  }
		public void setControlAcc(double controlAcc) {
		    super.setControlAcc(controlAcc);
		  }
		public int getN() {
		    return super.getN();
		  }
		public void setN(int n) {
		    super.setN(n);
		  }
		public double getTestAcc() {
		    return super.getTestAcc();
		  }
		public void setTestAcc(double testAcc) {
		    super.setTestAcc(testAcc);
		  }
		public double eval(Complex z) {
		
		    schottky.abelianDifferentialOf1stKind(controlValue, z, n, controlAcc);
		    schottky.abelianDifferentialOf1stKind(testValue, z, n, acc);
		
		    return controlValue.dist(testValue);
		  }
	}

	/**
   * Creates function which evaluates numerical the absolute error
   * we make when we use the approxiamation with a prescribed
   * error given by testAcc. This is simply done by measering the distance
   * of the approximation, with prescribed testAcc, to a second
   * approximation, with significant higher accuracy, the controlAcc
   * which is taken
   * is used for the approximation of differential of first kind.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @param n index of differential
   * @return function countin the number of elment in section.
   */
  public static SchottkyAnalysisFunction
      createNumericalErrorFunctionFor1stKind
      (Schottky schottky, double testAcc, double controlAcc, int n) {

    return new NumErr1stKind(schottky, testAcc, controlAcc, n);
  }

  /**
    * Creates function which evaluates numerical the absolute error
    * we make when we use the approxiamation with a prescribed
    * error given by testAcc, but this time we refer to the absolute series.
    * This is simply done by measering the distance
    * of the approximation, with prescribed testAcc, to a second
    * approximation, with significant higher accuracy, the controlAcc
    * which is taken
    * is used for the approximation of differential of first kind.
    * @param schottky instance of analysis
    * @param testAcc accuracy of approximation
    * @param n index of differential
    * @return function countin the number of elment in section.
    */
   public static SchottkyAnalysisFunction
       createNumericalErrorOfAbsFunctionFor1stKind
       (Schottky schottky, double testAcc, double controlAcc, int n) {

     return new NumericalErrorOfAbsFunctionFor1stKind(schottky, testAcc, controlAcc, n);
   }

    /**
     * Creates function counting the number of elements in subset which
     * is used for the approximation of the differential of first kind.
     * @param schottky instance of analysis
     * @param testAcc accuracy of approximation
     * @param n index of differential
     * @return function counting the number of elment in used subset.
   */
  public static SchottkyAnalysisFunction
      createNumOfUsedTermsFunctionFor1stKind
      (Schottky schottky, double testAcc, int n) {

    return new NumOfUsedTermsFunctionFor1stKind(schottky, testAcc, n);
  }


  /**
   * Creates function which delivers the smallest word length of all
   * element in section of the the subset which
   * is used for the approximation of differential of first kind.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @param n index of differential
   * @return function delivering smallest word length of section
   */
  public static SchottkyAnalysisFunction
      createMinWordLengthFunctionFor1stKind
      (Schottky schottky, double testAcc, int n) {

    return new MinWordLengthFunctionFor1stKind(schottky, testAcc, n);
  }

  /**
   * Creates function which checks the analysis object for
   * differential of the 1st kind.
   * This is done by comparing the result of the regulare
   * routine with the one for the analysis. This values must
   * coinside, otherwise something is wrong.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @param n index of differential
   * @return function delivering smallest word length of section
   */
  public static SchottkyAnalysisFunction
      createAnalysisCheckFunctionFor1stKind
      (Schottky schottky, double testAcc, int n) {

    return new AnalysisCheckFunctionFor1stKind(schottky, testAcc, n);
  }

  /**
   * Creates function which delivers the largest word length of all
   * element in section of the the subset which
   * is used for the approximation of differential of first kind.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @param n index of differential
   * @return function delivering largestest word length of section
   */
  public static SchottkyAnalysisFunction
      createMaxWordLengthFunctionFor1stKind
      (Schottky schottky, double testAcc, int n) {

    return new MaxWordLengthFunctionFor1stKind(schottky, testAcc, n);
  }

  /**
   * Creates function which delivers the absolute series
   * related to the differential of first kind.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @param n index of differential
   * @return function delivering the absolute series
   */
  public static SchottkyAnalysisFunction
      createAbsSeriesFunctionFor1stKind
      (Schottky schottky, double testAcc, int n) {

    return new AbsSeriesFunctionFor1stKind(schottky, testAcc, n);
  }

    /**
     * Creates function which delivers the theoretical error
     * of the approximation of differential of first kind.
     * @param schottky instance of analysis
     * @param testAcc accuracy of approximation
     * @param n index of differential
     * @return function countin the number of elment in section.
     */
    public static SchottkyAnalysisFunction
        createErrorFunctionFor1stKind
        (Schottky schottky, double testAcc, int n) {

      return new ErrorFunctionFor1stKind(schottky, testAcc, n);
    }


      /**
       * Creates function which delivers the waste caused by the
       * evaluation algorithm.
       * @param schottky instance of analysis
       * @param testAcc accuracy of approximation
       * @param n index of differential
       * @return function countin the number of elment in section.
       */
      public static SchottkyAnalysisFunction
          createWasteFunctionFor1stKind
          (Schottky schottky, double testAcc, int n) {

        return new WasteFunctionFor1stKind(schottky, testAcc, n);
      }

  /**
   * Creates function evaluating the absolute value of the
   * abelian differential of fist kind.
   * @param schottky instance of analysis
   * @param testAcc accuracy of approximation
   * @param n index of differential
   * @return function countin the number of elment in section.
   */
  public static SchottkyAnalysisFunction
      createAbsValueFunctionFor1stKind
      (Schottky schottky, double testAcc, int n) {

    return new AbsValueFunctionFor1stKind(schottky, testAcc, n);
  }


    /**
     * Creates function evaluating the argument value of the
     * abelian differential of fist kind.
     * @param schottky instance of analysis
     * @param testAcc accuracy of approximation
     * @param n index of differential
     * @return function countin the number of elment in section.
     */
    public static SchottkyAnalysisFunction
        createArgValueFunctionFor1stKind
        (Schottky schottky, double testAcc, int n) {

      return new ArgValueFunctionFor1stKind(schottky, testAcc, n);
    }

    /**
     * Creates function measering bench marks for provided function.
     * The bench mark is computed by evaluating 110 time the function
     * at given position. This is done in 10-packs. The one with the longest
     * duration of the elven 10-packs is excluded. Thus the bench mark is
     * finnally averaged over 100 evaluations.
     * @param function we want to banch mark.
     * @return bench mark function for analysis function.
     */
    public static SchottkyAnalysisFunction
        createBenchMarkFunction(SchottkyAnalysisFunction function) {

      return new BenchMarkFunction0(function);
    }

}
