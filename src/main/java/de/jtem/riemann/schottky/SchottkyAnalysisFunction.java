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
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;
import de.jtem.numericalMethods.calculus.minimizing.Info;
import de.jtem.numericalMethods.calculus.minimizing.Powell;

public abstract class SchottkyAnalysisFunction
    implements Serializable {

  public static final long serialVersionUID = 1L;

  Schottky schottky;

  double acc;
  double controlAcc;

  int n = 0;

  SchottkyAnalysisFunction(Schottky schottky, double testAcc) {
    this(schottky, testAcc, testAcc / 100, 0);
  }

  SchottkyAnalysisFunction(Schottky schottky, double testAcc, int n) {
    this(schottky, testAcc, testAcc / 100, n);
  }

  SchottkyAnalysisFunction(Schottky schottky,
                           double testAcc, double controlAcc) {
    this(schottky, testAcc, controlAcc, 0);
  }

  SchottkyAnalysisFunction(Schottky schottky,
                           double testAcc, double controlAcc, int n) {

    if (testAcc < controlAcc) {
      throw new IllegalArgumentException("test accuracy exeeds control");
    }

    this.schottky = schottky;
    this.acc = testAcc;
    this.controlAcc = controlAcc;
    this.n = n;
  }

  /**
   * returns instance which is to be analysed.
   */
  public Schottky getSchottky() {
    return schottky;
  }

  /**
   * returns index of function, e.g. index of differential.
   */
  int getN() {
    return n;
  }

  /**
   * sets index of function, e.g. index of differential.
   * @param n index of function, e.g. index of differential
   */
  void setN(int n) {
    this.n = n;
  }

  /**
   * returns accuracy of setup which is to analyse.
   */
  double getTestAcc() {
    return acc;
  }

  /**
   * sets accuracy of setup which is to analyse.
   * @param testAcc accuracy of setup
   */
  void setTestAcc(double testAcc) {
    this.acc = testAcc;
  }

  /**
   * returns accuracy of "precise" values which is used for reference.
   */
  double getControlAcc() {
    return controlAcc;
  }

  /**
   * sets accuracy of "precise" values which is used for reference.
   * This should atleast two digits preciser then testAcc.
       * @param controlAcc accuracy of "precise" values which is used for reference.
   */
  void setControlAcc(double controlAcc) {
    this.controlAcc = controlAcc;
  }

  final Complex testValue = new Complex();
  final Complex controlValue = new Complex();

  /**
   * Evaluates function at z.
   * @param z
   * @return function value
   */
  public abstract double eval(Complex z);

  /**
   * Evaluates function at all points provided thru xy and
   * stores result in valueAtXY.
   * @param xy point list: x0, y0, x1, y1, ....
   * @param valueAtXY value list: f( x0+iy0 ), f( x1+iy1 ), ...
   */
  public void eval(double[] xy, double[] valueAtXY) {
    Complex z = new Complex();

    for (int i = 0, j = 0; i < xy.length; j++) {
      z.assign(xy[i++], xy[i++]);
      valueAtXY[j] = eval(z);
    }
  }

  /**
   * Evaluates function at all points provided thru xy and
   * returns result in value list.
   * @param xy point list: x0, y0, x1, y1, ....
   * @return value list: f( x0+iy0 ), f( x1+iy1 ), ...
   */
  public double[] eval(double[] xy) {
    double[] valueAtXY = new double[xy.length / 2];
    eval(xy, valueAtXY);
    return valueAtXY;
  }

  public abstract static class ForProperty
      extends SchottkyAnalysisFunction {

    ForProperty(Schottky schottky, double testAcc) {
      super(schottky, testAcc, 0);
    }

    public double getTestAcc() {
      return super.getTestAcc();
    }

    public void setTestAcc(double testAcc) {
      super.setTestAcc(testAcc);
    }

  }

  public abstract static class ForIndexedProperty
      extends SchottkyAnalysisFunction.ForProperty {

    ForIndexedProperty(Schottky schottky, double testAcc, int n) {
      super(schottky, testAcc);

      super.n = n;
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

  }

  /**
   * Blue print for function computing bench marks for an other.
   */
  static abstract class BenchMarkFunction
      extends SchottkyAnalysisFunction {

    SchottkyAnalysisFunction function;

    BenchMarkFunction(SchottkyAnalysisFunction function) {
      super(function.schottky, function.acc, function.n);
      this.function = function;
    }
  }

  abstract class ExtremizerFunctional
      implements RealFunctionOfSeveralVariables {

    Complex z = new Complex();

    public int getNumberOfVariables() {
      return 2;
    }
  }

  class MinimizerFunctional
      extends ExtremizerFunctional {

    public double eval(double p[]) {
      z.assign(p[0], p[1]);

      if (schottky.isInFundamentalDomain(z, 1e-12)) {
        return SchottkyAnalysisFunction.this.eval(z);
      }

      return Double.MAX_VALUE;
    }
  }

  class MaximizerFunctional
      extends ExtremizerFunctional {

    public double eval(double p[]) {
      z.assign(p[0], p[1]);

      if (schottky.isInFundamentalDomain(z, 1e-12)) {
        return 1 / SchottkyAnalysisFunction.this.eval(z);
      }

      return Double.MAX_VALUE;
    }

  }

  double extremize(Complex z, double ftol,
                   ExtremizerFunctional f, int n, boolean debug) {

    Complex originalZ = new Complex(z);

    double[] value = new double[] {
        z.re, z.im};

    Powell.search(value, ftol, n, f, new Info(debug));

    z.assign(value[0], value[1]);

    return eval(z);
  }

  /**
   * Minimizes function using powells algorith starting at z.
   * @param z start postion of minimization and minimum on output.
   * @return value of minimium
   */
  public double minimize(Complex z) {
    return minimize(z, 50, 0);
  }

  /**
   * Performs a minimizing step of powells algorithm starting at z.
   * @param z start position of minimization and minimum on output.
   * @return value of minimium
   */
  public double minimizeStep(Complex z) {
    return minimize(z, 1e-10, 1);
  }

  /**
   * Minimizes function using powells algorithm starting at z.
   * @param z start postion of minimization and minimum on output.
   * @param ftol tolerance of powells algorithm
   * @param n maximal number of steps in powells algorithm
   * @return value of minimium
   */
  public double minimize(Complex z, double ftol, int n) {
    return minimize(z, ftol, n, false);
  }

  /**
   * Minimizes function using powells algorithm starting at z.
   * @param z start postion of minimization and minimum on output.
   * @param ftol tolerance of powells algorithm
   * @param n maximal number of steps in powells algorithm
   * @param toggle for verbose
   * @return value of minimium
   */
  public double minimize(Complex z, double ftol,
                         int n, boolean debug) {

    return extremize(z, ftol, new MinimizerFunctional(), n, debug);
  }

  /**
   * Maximizes function using powells algorith starting at z.
   * @param z start postion of maximization and maximum on output.
   * @return value of maximium
   */
  public double maximize(Complex z) {
    return maximize(z, 50, 0);
  }

  /**
   * Performs a maximizing step of powells algorithm starting at z.
   * @param z start position of maximization and maximum on output.
   * @return value of maximium
   */
  public double maximizeStep(Complex z) {
    return maximize(z, 1e-10, 1);
  }

  /**
   * Maximizes function using powells algorithm starting at z.
   * @param z start postion of maximization and maximum on output.
   * @param ftol tolerance of powells algorithm
   * @param n maximal number of steps in powells algorithm
   * @return value of maximium
   */
  public double maximize(Complex z, double ftol, int n) {
    return maximize(z, ftol, n, false);
  }

  /**
   * Maximizes function using powells algorithm starting at z.
   * @param z start postion of maximization and maximum on output.
   * @param ftol tolerance of powells algorithm
   * @param n maximal number of steps in powells algorithm
   * @param toggle for verbose
   * @return value of maximium
   */
  public double maximize(Complex z, double ftol,
                         int n, boolean debug) {

    return extremize(z, ftol, new MaximizerFunctional(), n, debug);
  }

}