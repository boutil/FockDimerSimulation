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

package de.jtem.riemann.theta;

import java.io.Serializable;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.blas.RealMatrix;
import de.jtem.blas.RealVector;
import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.algebra.linear.MatrixOperations;
import de.jtem.numericalMethods.algebra.linear.decompose.Cholesky;
import de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;
import de.jtem.numericalMethods.geometry.latticeReduction.LLL;

/**
 * Riemann theta function
 * <p>
 * The de.jtem.riemann theta function is defined by its mulit-dimensional Fourier series
 * <p align=center>
 *   <code>
 *     &theta;(z|B) = &sum;<sub>n<small>&isin;</small>Z<sup>g</sup></sub> exp( &frac12;(Bn,n)+(z,n) )
 *   </code>,
 * </p>
 * where
 * <code>z<small>&isin;</small>C<sup>g</sup>,B<small>&isin;</small>C<sup>g<small>x</small>g</sup></code>.
 * <code>B</code> is symmetric and has strictly negative definite real part.
 * <p>
 * The exponential growth in <code>z</code> of the Riemann theta function can be factored
 * out and the remaining oscillatory part can be approximated with a prescribed error
 * tolerance. Thus we have:
 * <p align=center>
 *   <code>
 *     &theta;(z|B) = exp( factor(z|B) ) &middot; thetaSum(z|B)
 *   </code>,
 * </p>
 * A detailed description can be found in the paper <em>Computing Riemann Theta Functions</em>,
 * which provides the theoretical background of this implementation.
 * <p>
 * Separating the exponential factor is crucial when dealing with theta functions. In practice,
 * they often appear as ratios, in which case the exponential growth can be canceled.
 * For example,
 * let <code>a,b,c,d <small>&isin;</small> C<sup>g</sup></code> with
 * <code>a+b = c+d mod 2&pi;iZ<sup>g</sup></code>. Then
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td>              </td> <td> &theta;(z+a|B)&middot;&theta;(z+b|B) </td> </tr>
 * <tr> <td>f(z) = &nbsp; </td> <td> <hr noshade size=1>                  </td> </tr>
 * <tr> <td>              </td> <td> &theta;(z+c|B)&middot;&theta;(z+d|B) </td> </tr>
 * </table>
 * </code>
 * <p>
 * defines an Abelian function.
 * The following code fragment is taken from {@link AbelianFunction},
 * it shows how you evaluate <code>f(z)</code> using <code>Theta</code>:
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td> Theta theat = new Theta( B, tol );                          </td></tr>
 * <tr> <td>                                                             </td></tr>
 * <tr> <td> theta.theta( z.plus(a), factorOfZPlusA, thetaSumOfZPlusA ); </td></tr>
 * <tr> <td> theta.theta( z.plus(b), factorOfZPlusB, thetaSumOfZPlusB ); </td></tr>
 * <tr> <td> theta.theta( z.plus(c), factorOfZPlusC, thetaSumOfZPlusC ); </td></tr>
 * <tr> <td> theta.theta( z.plus(d), factorOfZPlusD, thetaSumOfZPlusD ); </td></tr>
 * <tr> <td>                                                             </td></tr>
 * <tr> <td> factor.assign(      factorOfZPlusA );                       </td></tr>
 * <tr> <td> factor.assignPlus(  factorOfZPlusB );                       </td></tr>
 * <tr> <td> factor.assignMinus( factorOfZPlusC );                       </td></tr>
 * <tr> <td> factor.assignMinus( factorOfZPlusD );                       </td></tr>
 * <tr> <td>                                                             </td></tr>
 * <tr> <td> fOfZ.assignExp( factor );                                   </td></tr>
 * <tr> <td> fOfZ.assignTimes(  thetaSumOfZPlusA );                      </td></tr>
 * <tr> <td> fOfZ.assignTimes(  thetaSumOfZPlusB );                      </td></tr>
 * <tr> <td> fOfZ.assignDivide( thetaSumOfZPlusC );                      </td></tr>
 * <tr> <td> fOfZ.assignDivide( thetaSumOfZPlusD );                      </td></tr>
 * </table>
 * </code>
 * <p>
 * Evaluating the theta function may easily lead to overflows or numerical
 * instabilities, which should be avoided. ( The code above avoids the creation of
 * temporarily used objects, which can be very beneficial for the performance. )
 * The first line of the code creates a theta function for a period matrix
 * <code>B</code> and a tolerance of <code>tol</code> for the <code>thetaSum(z|B)</code>.
 * This results in an absolute error of
 * <p align=center>
 * <code>
 * <table>
 * <tr><td>thetaSumZPlusA+thetaSumZPlusB        </td> <td>                                   </td></tr>
 * <tr><td><hr noshade size=1>                  </td> <td> &middot; exp(factor) &middot; tol </td></tr>
 * <tr><td>thetaSumZPlusC&middot;thetaSumZPlusD </td> <td>                                   </td></tr>
 * </table>
 * </code>
 * <p>
 * for <code>fOfZ</code>. This error has be  be controled if you want to be really shure what
 * is going on. In practice it suffies to avoid the poles of <code>f</code>. This is true because
 * the estimates, which lead to the approximation of the oscillating part of the Riemann
 * theta function, are rather weak, and therefore the error for the actual values often falls low
 * of the tolerance. This is also the reason why we have never experienced
 * inaccuracies when evaluating first and second order derivatives of the Riemann
 * theta function. Hence, we have not incorporated error estimates for these.
 * <p>
 * To create a theta function object you must provide the period matrix and you may provide
 * an error tolerance. Additionally, you may switch off features, which 
 * are used by default.
 * <p>
 * The first is about the fill factor error. The fill factor error is a heuristic to improve
 * the accuracy of the error prediction. As mentioned earlier the estimates which lead to
 * the approximation error are rather coarse, which is improved by this heuristic for
 * the price of loosing the certainty of the 100% error. The smaller
 * the fill factor, the bigger its influence on the approximation, which
 * anyway is usually irrelevant for small genus (<code>g<4</code>).
 * <p>
 * The second feature is Siegel`s reduction. <code>{@link SiegelReduction}</code> uses the modular
 * transformation property of the Riemann theta function and can reduce the cost of
 * computing theta functions drasticly. Sometimes you already know that your period
 * matrix cannot be reduced. In that case you would switch of this feature to save the effort
 * of checking for a possible reduction.
 * The <code>{@link ModularPropertySupport}</code> is implemented without controling
 * a nasty 8<sup>th</sup>-root of unity,
 * which leads to an ambiguity. It is, however, usually neglectable because one only deals
 * with ratios of theta functions. In case you are interested in theta functions themselves
 * you still can take advantage of Siegel`s reduction: just evaluate theta at zero
 * (or any other point) with and
 * without Siegel`s reduction and determine that 8<sup>th</sup>-root of unity.
 * <p>
 * Finally you can decide whehter you want to use pointwise or uniform approximation.
 * The uniform approximation is the default. Allthough we experienced that for genus g>3
 * the pointwise is superior. If you experience memory problems, because of the number
 * of terms in the theta sum, you should toggle off the uniform approximation and use the
 * the pointwise.
 * <p>
 * Further reading:
 * <table>
 * <tr><td>-</td><td>B.Deconinck, M.Heil, A.I.Bobenko, M.van Hoeij, and M.Schmies.
 * <em>Computing Riemann Theta Functions</em>. Is to appear. </td></tr>
 * <tr><td>-</td><td>M.Heil.
 * <em>Numerical Tools for the study of finite gap solutions of integrable systems.</em>
 * PhD thesis, Technische Universit&auml;t Berlin, 1995. </td></tr>
 * <tr><td>-</td><td>E.D.Belokolos, A.I.Bobenko, V.Z.Enolskii, A.R.Its, and V.B.Matveev.
 * <em>Algebro-geometric appoach to nonlinear integrable probems.</em>
 * Springer-Verlag Berlin, 1994. </td></tr>
 * <tr><td>-</td><td>C.L.Siegel.
 * <em>Topics in complex function theory. Vol.III.</em>
 * John Wiley & Sons, Inc., New York, 1989</td></tr>
 * </table>
 * @see ThetaWithChar
 * @see TransformPropertySupport
 * @see ModularPropertySupport
 * @see SiegelReduction
 * @see AbelianFunction
 * @author Markus Schmies
 */

public class Theta extends AbstractTheta implements Serializable, Cloneable {

	private static final long serialVersionUID = 1L;

	private static final Complex PI2I = new Complex(0, 2 * Math.PI);

	double tol = 1E-7;

	int dim;

	ComplexMatrix periodMatrix = new ComplexMatrix();

	TransformPropertySupport transform;
	ModularPropertySupport modular;
	SiegelReduction siegel;

	boolean performSiegelReduction = true;
	boolean useFillFactorError = true;
	boolean uniformApproximation = true;
	
	long lastChangeOfPeriodMatrix = System.currentTimeMillis();

	//LatticePointsForUniformApproximation latticePointsforUniformApproximation;
	Object latticePointsforApproximation;
	
	RealMatrix latticePoints;  // only used in case of uniform approximation
	int numOfLatticePoints;  // only used in case of uniform approximation

	ComplexMatrix B;

	RealMatrix reB = new RealMatrix();
	RealMatrix imB = new RealMatrix();
	RealMatrix T = new RealMatrix();

	ComplexVector expOfHalfBnn;

	double radius = 1;
	double lSLV; // length of shortest lattice vector
	double fillFactor;
	double detOfLattice;
	double lSLV2PowOfDim; // = Math.pow( lSLV, dim )
	double two2PowOfDim; // = Math.pow( 2,    dim )
	double gammaOfHalfDim; // = Gamma.of( dim / 2.0 )
	double sqrtOfPI2PowOfDim; // = Math.pow( Math.sqrt( Math.PI ), dim )

	int highestOrder = 0;

	boolean modularIsId;

	ErrorRadiusSolverFunction ersf = new ErrorRadiusSolverFunction();

	/* the following latticePointss are used as temporary variables in
	   varias methods which makes the whole class not thread save!!! */

	private final Complex tmp = new Complex();

	private final Complex continousCorrection = new Complex();
	private final Complex exponent = new Complex();
	private final Complex term = new Complex();
	private final Complex termX = new Complex();
	private final Complex termY = new Complex();
	private final Complex termXY = new Complex();
	private final Complex nX = new Complex();
	private final Complex nY = new Complex();
	private final Complex nZ = new Complex();
	private final Complex expOfNZ = new Complex();
	private final Complex invOfExpOfNZ = new Complex();

	private final Complex MX = new Complex();
	private final Complex MY = new Complex();

	private final Complex factor = new Complex();
	private final Complex dFactorByDX = new Complex();
	private final Complex dFactorByDY = new Complex();
	private final Complex ddFactorByDXDY = new Complex();
	private final ComplexVector dTByDX = new ComplexVector();
	private final ComplexVector dTByDY = new ComplexVector();
	private final ComplexVector ddTByDXDY = new ComplexVector();

	private final ComplexVector tmpVector = new ComplexVector();

	/**
	 * Creates a Riemann theta function with prescribed <code>periodMatrix</code>.
	 * The error tolerance is by default <code>e-7</code>, the fill factor error
	 * and Siegel`s reduction is used.
	 * @param periodMatrix symmetric complex matrix with negative definite real part
	 */
	public Theta(final ComplexMatrix periodMatrix) {
		this(periodMatrix, 1e-7, true, true);
	}

	/**
	 * Creates a Riemann theta function with prescribed <code>periodMatrix</code>
	 * and error tolerance <code>tol</code>.
	 * Both, fill factor error and Siegel`s reduction is used.
	 * @param periodMatrix symmetric complex matrix with negative definite real part
	 * @param tol positive number
	 */
	public Theta(final ComplexMatrix periodMatrix, final double tol) {
		this(periodMatrix, tol, true, true);
	}

	/**
	 * Creates a Riemann theta function with prescribed <code>periodMatrix</code>
	 * and error tolerance <code>tol</code>.
	 * It further enables the user to switch off the Siegel`s Reduction.
	 * @param periodMatrix symmetric complex matrix with negative definite real part
	 * @param tol positive number
	 * @param performSiegelReduction  controles whether or not Siegel`s reduction is performed
	 */
	public Theta(final ComplexMatrix periodMatrix, final double tol, final boolean performSiegelReduction) {
		this(periodMatrix, tol, performSiegelReduction, true);
	}

	/**
	 * Creates a Riemann theta function
	 * with prescribed <code>periodMatrix</code>
	 * and error tolerance <code>tol</code>.
	 * It further enables the user to switch off the Siegel`s Reduction and the
	 * heuristic of the fill factor error.
	 * @param periodMatrix symmetric complex matrix with negative definite real part
	 * @param tol positive number
	 * @param performSiegelReduction  controles whether or not Siegel`s reduction is performed
	 * @param useFillFactorError controles whether the fill factor or the 100% error is used.
	 */
	public Theta(final ComplexMatrix periodMatrix, final double tol, final boolean performSiegelReduction, final boolean useFillFactorError) {

		this.performSiegelReduction = performSiegelReduction;

		this.tol = tol;

		setPeriodMatrix(periodMatrix);
	}

	/**
	 * Creates a Riemann theta function
	 * with prescribed <code>periodMatrix</code>
	 * and error tolerance <code>tol</code>.
	 * It further enables the user to switch off the Siegel`s Reduction and the
	 * heuristic of the fill factor error.
	 * @param periodMatrix symmetric complex matrix with negative definite real part
	 * @param tol positive number
	 * @param performSiegelReduction  controles whether or not Siegel`s reduction is performed.
	 * @param useFillFactorError controles whether the fill factor or the 100% error is used.
	 * @param uniformApproximation controles whether to use uniform or pointwise approximation.
	 */
	public Theta(final ComplexMatrix periodMatrix, final double tol, final boolean performSiegelReduction, final boolean useFillFactorError, final boolean uniformApproximation ) {

		this.performSiegelReduction = performSiegelReduction;
		this.useFillFactorError = useFillFactorError;
		this.uniformApproximation = uniformApproximation;
		
		this.tol = tol;

		setPeriodMatrix(periodMatrix);
	}
	
	/**
	 * Returns the fill factor.
	 */
	public final double getFillFactor() {
		return fillFactor;
	}

	/**
	 * Returns whether the heuristic of the fill factor error is used or not.
	 */
	public final boolean getUseFillFactorError() {
		return useFillFactorError;
	}

	/**
	 * Switches the heuristic fill factor error on or off.
	 * @param useFillFactorError controles whether the fill factor or the 100% error is used.
	 */
	public final void setUseFillFactorError(final boolean useFillFactorError) {
		if (useFillFactorError == this.useFillFactorError)
			return;

		this.useFillFactorError = useFillFactorError;

		computeRadius();
		computeLatticePoints();
	}

	/**
	 * Returns whether  Siegel's reduction algorithm is peroformed or not.
	 */
	public final boolean isSiegelReductionPerformed() {
		return performSiegelReduction;
	}

	/**
	 * Switches between pointwise and uniform approximation.
	 * @param performSiegelReduction controles whether  Siegel's reduction algorithm is peroformed or not.
	 */
	public final void setSiegelReductionPerformed(final boolean performSiegelReduction) {
		if (performSiegelReduction == this.performSiegelReduction)
			return;

		this.performSiegelReduction = performSiegelReduction;

		compute();
	}
	
	/**
	 * Returns whether pointwise or uniform approximation is used.
	 */
	public final boolean isUniformApproximation() {
		return uniformApproximation;
	}

	/**
	 * Switches between pointwise and uniform approximation.
	 * @param uniformApproximation controles whether pointwise or uniform approximation is used.
	 */
	public final void setUniformApproximation(final boolean uniformApproximation) {
		if (uniformApproximation == this.uniformApproximation)
			return;

		this.uniformApproximation = uniformApproximation;

		this.latticePointsforApproximation = null;
		
		computeLatticePoints();
	}
	
	/**
	 * Returns accuracy for the oscillatory part of the theta function.
	 */
	public final double getAccuracy() {
		return tol;
	}

	/**
	 * Sets the accuracy for the oscillatory part of the theta function.
	 * @param accuracy positive number
	 */
	public final void setAccuracy(final double accuracy) {
		if (this.tol == accuracy)
			return;

		this.tol = accuracy;

		computeRadius();
		computeLatticePoints();
	}
	
	/**
	 *  @deprecated use {@link #getAccuracy()} 
	 * Returns the error tolerance for the oscillatory part of the theta function.
	 */
	public final double getTol() {
		return getAccuracy();
	}

	/**
	 * @deprecated use {@link #setAccuracy(double)} 
	 * Sets the error tolerance for the oscillatory part of the theta function.
	 * @param tol positive number
	 */
	public final void setTol(final double tol) {
		setAccuracy( tol );
	}

	/** Returns the dimension of the complex vector space on which the Riemann theta
	functions operates.
	*/
	public final int getDim() {
		return dim;
	}

	public ModularPropertySupport getModularPropertySupport() {
		return modular;
	}
	
	void setDim(final int v) {

		if (v != dim) {

			dim = v;

			transform = new TransformPropertySupport(dim);

			siegel = new SiegelReduction(dim);
			modular = new ModularPropertySupport(dim);

			gammaOfHalfDim = Gamma.gamma(dim / 2.0);
			two2PowOfDim = Math.pow(2, dim);
			sqrtOfPI2PowOfDim = Math.pow(Math.sqrt(Math.PI), dim);
		}
	}

	/**
	 * Returns number of lattice points which are used for the uniform approximation
	 * of the oscillatory part of the Riemann theta function.
	 */
	public final int getNumOfLatticePoints() {
		return numOfLatticePoints;
	}

	/**
	 * Returns period matrix on which the Riemann theta function operates.
	 */
	public final ComplexMatrix getPeriodMatrix() {
		return new ComplexMatrix(periodMatrix);
	}

	public final ComplexMatrix getB() {
		return new ComplexMatrix(B);
	}

	/**
	 * Sets periodMatrix on which the Riemann theta function operates.
	 * @param periodMatrix symmetric complex matrix with negative definite real part
	 */
	public final void setPeriodMatrix(final ComplexMatrix periodMatrix) {
		if (this.periodMatrix != null && periodMatrix.equals(this.periodMatrix))
			return;

		if (!periodMatrix.isSquared())
			throw new IllegalArgumentException("matrix is not squared");
		if (!periodMatrix.isSymmetric())
			throw new IllegalArgumentException("matrix is not symmetric");

		setDim(periodMatrix.getNumRows());

		if (this.periodMatrix == null)
			this.periodMatrix = new ComplexMatrix(periodMatrix);
		else
			this.periodMatrix.assign(periodMatrix);

		lastChangeOfPeriodMatrix = System.currentTimeMillis();

		compute();
	}

	final void compute() {

		if (performSiegelReduction && dim >= 1) {

			siegel.setPeriodMatrix(periodMatrix);

			modular.setModularTransformation(siegel.modular);

			B = siegel.tPM;

			T.assignTranspose(siegel.tL); // -re(B) = tL * T;

			// the lattice:  sqrt(1/2) T * Z^g
			detOfLattice = T.get(0, 0);
			for (int i = 1; i < dim; i++)
				detOfLattice *= T.get(i, i);
			detOfLattice /= Math.sqrt(two2PowOfDim);

			lSLV = T.get(0, 0) / Math.sqrt(2);
		} else {

			modular.assignId();

			B = periodMatrix;

			B.getRe(T);

			T.assignTimes(-1);

			Cholesky.decompose(T.re, T.re);

			T.assignTranspose();

			// the lattice:  sqrt(1/2) T * Z^g
			detOfLattice = T.get(0, 0);
			for (int i = 1; i < dim; i++)
				detOfLattice *= T.get(i, i);
			detOfLattice /= Math.sqrt(two2PowOfDim);

			T.assignTranspose();

			LLL.reduce(T.re, (int[][]) null);

			lSLV = 0;
			for (int i = 0; i < dim; i++)
				lSLV += T.get(0, i) * T.get(0, i);
			lSLV = Math.sqrt(lSLV /= 2);
		}

		lSLV2PowOfDim = Math.pow(lSLV, dim);

		fillFactor = 2 * sqrtOfPI2PowOfDim * lSLV2PowOfDim / two2PowOfDim / dim / gammaOfHalfDim / detOfLattice;

		B.getRe(reB);
		B.getIm(imB);

		modularIsId = modular.isId();

		if (!modularIsId) {
			modular.setPeriodMatrix(periodMatrix);
		}

		transform.setPeriodMatrix(B);

		computeRadius();
		computeLatticePoints();
		
	}

	/**
	 * Returns radius of ellipsoid which is used to determine the lattice
	 * points which are used for the uniform approximation
	 * of the oscillatory part of the Riemann theta function.
	 */
	public final double getRadius() {
		return radius;
	}

	
	void computeLatticePoints() {

		if( uniformApproximation ) {
			computeLatticePointsForUniformApproximation();
		} else {
			computeLatticePointsForPointwiseApproximation();
		}
	
	}
	
	private void computeLatticePointsForPointwiseApproximation() {
		
		LatticePointsInEllipsoidIterator iterator = 
			(LatticePointsInEllipsoidIterator)this.latticePointsforApproximation;
		
		if( iterator == null ) {
			this.latticePointsforApproximation =
				iterator = new LatticePointsInEllipsoidIterator( reB.times(-0.5) );
		} else {
				iterator.setB(reB.times(-0.5));
		}
	}

	private void computeLatticePointsForUniformApproximation() {
		
		LatticePointsForUniformApproximation 
			latticePointsforUniformApproximation = (LatticePointsForUniformApproximation)this.latticePointsforApproximation;
		// the lattice:  sqrt(1/2) T * Z^g
		if (latticePointsforUniformApproximation == null) {
			//latticePointsforUniformApproximation = new LatticePointsForUniformApproximation(reB.times(-0.5), true);
			this.latticePointsforApproximation =
				latticePointsforUniformApproximation = 
					new LatticePointsForUniformApproximation(reB.times(-0.5));
		} else {
			latticePointsforUniformApproximation.setB( reB.times(-0.5)  );
		}
		
		latticePointsforUniformApproximation.setRadius(radius);

		latticePointsforUniformApproximation.update();

		latticePoints = latticePointsforUniformApproximation.latticePoints;

		numOfLatticePoints = latticePointsforUniformApproximation.numOfLatticePoints / 2 + 1;
		
		// the thetaSum procedures expect the zero at the beginning
		double [] tmp = latticePoints.re[numOfLatticePoints-1];
		latticePoints.re[numOfLatticePoints-1] =  latticePoints.re[0];;
		latticePoints.re[0] = tmp;
	
		if (expOfHalfBnn == null || expOfHalfBnn.size() < numOfLatticePoints) {
			expOfHalfBnn = new ComplexVector(numOfLatticePoints);
		}

		final double[][] BRe = B.re;
		final double[][] BIm = B.im;

		for (int i = 0; i < numOfLatticePoints; i++) {

			final double[] nRe = latticePoints.re[i];

			double re = 0;
			double im = 0;

			for (int j = 0; j < dim; j++) {

				final double[] rowBRe = BRe[j];
				final double[] rowBIm = BIm[j];

				re += rowBRe[j] * nRe[j] * nRe[j] / 2;
				im += rowBIm[j] * nRe[j] * nRe[j] / 2;

				for (int k = j + 1; k < dim; k++) {

					re += rowBRe[k] * nRe[j] * nRe[k];
					im += rowBIm[k] * nRe[j] * nRe[k];
				}
			}

			exponent.assignExp(re, im);

			expOfHalfBnn.set(i, exponent);
		}
	}

	final void error0(double x, double[] value) {

		double f = dim * two2PowOfDim / lSLV2PowOfDim / 2.0;

		if (useFillFactorError) {
			f *= fillFactor;
		}

		value[0] = x;
		value[1] = f * Gamma.gamma(dim / 2.0, x, gammaOfHalfDim);
		value[2] = -1 * f * Math.exp(Math.log(x) * (dim / 2.0 - 1) - x);
	}

	final class ErrorRadiusSolverFunction implements NewtonRaphson.RealFunctionWithDerivative, Serializable {
		double error;

		int order = 0;

		double[] rootValues = new double[3];

		public void eval(double x, double[] f, int offsetF, double[] df, int offsetDF) {

			switch (order) {
				case 1 :
					break;
				case 2 :
					break;
				default :

					double factor = dim * two2PowOfDim / lSLV2PowOfDim / 2.0;

					if (useFillFactorError) {
						factor *= fillFactor;
					}

					f[offsetF] = factor * Gamma.gamma(dim / 2.0, x, gammaOfHalfDim);
					df[offsetDF] = -1 * factor * Math.exp(Math.log(x) * (dim / 2.0 - 1) - x);
			}

			f[offsetF] -= error;
		}

		double getRadius(double anError, int anOrder) {

			error = anError;
			order = anOrder;

			NewtonRaphson.search(this, dim / 2.0, rootValues);

			return lSLV / 2 + Math.sqrt(rootValues[0]);
		}
	}

	/**
	 * Returns shortest lattice vector of possibly Siegel reduced period matrix.
	 * @return shortest lattice vector
	 */
	public double getShortestLatticeVector() {
		return lSLV;
	}
	
	final void computeRadius() {

		double harmonicThreshRadius = (Math.sqrt(dim + 2 * highestOrder + Math.sqrt(dim * dim + 8 * highestOrder)) + lSLV) / 2;

		try {

			radius = Math.max(ersf.getRadius(tol, 0), harmonicThreshRadius);

		} catch (RuntimeException e) {

			radius = harmonicThreshRadius;
		}
	}

	private RealVector x = new RealVector();
	private RealVector y = new RealVector();
	
	private RealVector c = new RealVector();
	
	private RealVector Xn = new RealVector();
	private RealVector Yn = new RealVector();
	private RealVector n = new RealVector();
	
	void thetaSumPointwise(final ComplexVector Z, final Complex thetaSumZ) {
		
		LatticePointsInEllipsoidIterator iterator = 
			(LatticePointsInEllipsoidIterator)this.latticePointsforApproximation;
		
		thetaSumZ.assign(0);
		
		numOfLatticePoints=0;
//		
//		this.n.newSize(dim);
//		this.n.re = iterator.n; //hack
//		
//		x.assign( Z.re );
//		y.assign( Z.im );
//		
//		c.assignTimes( transform.reBInv, x); 
//		c.assignNeg();
		
		c.newSize(dim);
		MatrixOperations.times(transform.reBInv.re,Z.re,c.re);
		c.assignNeg();
		
		final double [] n = iterator.n;
		final double [] x = Z.re;
		final double [] y = Z.im;
		
		iterator.startIteration( radius, c.re );
		
		while( iterator.hasNext() ) {
			numOfLatticePoints++;			
			
			double nXn = 0;
			double nYn = 0; 
			
			double nx = 0;
			double ny = 0;
			
			for( int i=0; i<dim; i++ ) {
				double Xn = 0;
				double Yn = 0;
				
				final double [] rowX = reB.re[i];
				final double [] rowY = imB.re[i];
				
				for( int j=0; j<dim; j++ ) {
					Xn+=rowX[j] *n[j];
					Yn+=rowY[j] *n[j];
				}
				
				nXn += n[i]*Xn;
				nYn += n[i]*Yn;
				
				nx += n[i]*x[i];
				ny += n[i]*y[i];
			}
//			Xn.assignTimes( reB, this.n );
//			Yn.assignTimes( imB, this.n );
//			
//			final double _nXn = this.n.dot(Xn);
//			final double _nYn = this.n.dot(Yn);
//			
//			final double _nx = this.n.dot(this.x);
//			final double _ny = this.n.dot(this.y);
			
			term.assignExp( nXn/2 + nx, nYn/2 +ny );
			
			thetaSumZ.assignPlus( term );
		}
	}
	
	void thetaSum(final ComplexVector Z, final Complex thetaSumZ) {
		if( this.uniformApproximation ) {
			thetaSumUniform( Z, thetaSumZ );
		} else {
			thetaSumPointwise( Z, thetaSumZ );
		}
	}
	
	void thetaSumUniform(final ComplexVector Z, final Complex thetaSumZ) {

		final double[] zRe = Z.re, zIm = Z.im;

		final double[][] intLatticePointsRe = latticePoints.re;

		final double[] expOfHalfBnnRe = expOfHalfBnn.re;
		final double[] expOfHalfBnnIm = expOfHalfBnn.im;

		thetaSumZ.assign(1, 0);

		for (int i = 1; i < numOfLatticePoints; i++) {

			final double[] nRe = intLatticePointsRe[i];

			double nZRe = 0, nZIm = 0;

			for (int j = 0; j < dim; j++) {

				nZRe += zRe[j] * nRe[j];
				nZIm += zIm[j] * nRe[j];
			}

			expOfNZ.assignExp(nZRe, nZIm);
			invOfExpOfNZ.assignInvert(expOfNZ);

			// compute thetaSumZ
			term.assignPlus(expOfNZ, invOfExpOfNZ);
			term.assignTimes(expOfHalfBnn.re[i], expOfHalfBnn.im[i]);
			thetaSumZ.assignPlus(term);
		}
	}

	
	void dThetaSumPointwise(final ComplexVector Z, final ComplexVector X, final Complex thetaSumZ, final Complex thetaSumX) {
			
		LatticePointsInEllipsoidIterator iterator = 
			(LatticePointsInEllipsoidIterator)this.latticePointsforApproximation;
		
		thetaSumZ.assign(0);
		thetaSumX.assign(0);
		
		numOfLatticePoints=0;
		
		n.newSize(dim);
		n.re = iterator.n; //hack
		
		x.assign( Z.re );
		y.assign( Z.im );
		
		c.assignTimes( transform.reBInv, x); 
		c.assignNeg();
		
		iterator.startIteration( radius, c.re );
		
		while(iterator.hasNext()) {
			numOfLatticePoints++;			
			
			Xn.assignTimes( reB, n );
			Yn.assignTimes( imB, n );
			
			final double nXn = n.dot(Xn);
			final double nYn = n.dot(Yn);
			
			final double nx = n.dot(x);
			final double ny = n.dot(y);
			
			term.assignExp( nXn/2 + nx, nYn/2 +ny );
			
			ComplexVector.dotBilinear( n, X, nX);
			
			termX.assignTimes( term, nX);
			
			thetaSumZ.assignPlus( term );
			thetaSumX.assignPlus( termX );
	
		}
	}
	
	void dThetaSum(final ComplexVector Z, final ComplexVector X, final Complex thetaSumZ, final Complex thetaSumX) {
		if( this.uniformApproximation ) {
			dThetaSumUniform( Z, X, thetaSumZ, thetaSumX );
		} else {
			dThetaSumPointwise( Z, X, thetaSumZ,  thetaSumX );
		}
	}
	
	final void dThetaSumUniform(final ComplexVector Z, final ComplexVector X, final Complex thetaSumZ, final Complex thetaSumX) {

		final double[] zRe = Z.re, zIm = Z.im;
		final double[] uRe = X.re, uIm = X.im;

		final double[][] intLatticePointsRe = latticePoints.re;

		final double[] expOfHalfBnnRe = expOfHalfBnn.re, expOfHalfBnnIm = expOfHalfBnn.im;

		thetaSumZ.assign(1, 0);
		thetaSumX.assign(0, 0);

		for (int i = 1; i < numOfLatticePoints; i++) {

			final double[] nRe = intLatticePointsRe[i];

			double nZRe = 0, nZIm = 0;
			double nXRe = 0, nXIm = 0;

			for (int j = 0; j < dim; j++) {
				final double n = nRe[j];

				nZRe += zRe[j] * n;
				nZIm += zIm[j] * n;

				nXRe += uRe[j] * n;
				nXIm += uIm[j] * n;
			}

			expOfNZ.assignExp(nZRe, nZIm);
			invOfExpOfNZ.assignInvert(expOfNZ);

			// compute thetaSumZ
			term.assignPlus(expOfNZ, invOfExpOfNZ);
			term.assignTimes(expOfHalfBnn.re[i], expOfHalfBnn.im[i]);
			thetaSumZ.assignPlus(term);

			// compute thetaSumX
			term.assignMinus(expOfNZ, invOfExpOfNZ);
			term.assignTimes(expOfHalfBnn.re[i], expOfHalfBnn.im[i]);
			term.assignTimes(nXRe, nXIm);

			thetaSumX.assignPlus(term);
		}
	}

	
	void ddThetaSumPointwise( 
			final ComplexVector Z,
			final ComplexVector X,
			final ComplexVector Y,
			final Complex thetaSumZ,
			final Complex thetaSumX,
			final Complex thetaSumY,
			final Complex thetaSumXY) {
		
		LatticePointsInEllipsoidIterator iterator = 
			(LatticePointsInEllipsoidIterator)this.latticePointsforApproximation;
		
		thetaSumZ.assign(0);
		thetaSumX.assign(0);
		thetaSumY.assign(0);
		thetaSumXY.assign(0);
		
		numOfLatticePoints=0;
		
		n.newSize(dim);
		n.re = iterator.n; //hack
		
		x.assign( Z.re );
		y.assign( Z.im );
		
		c.assignTimes( transform.reBInv, x); 
		c.assignNeg();
		
		iterator.startIteration( radius, c.re );
		while(iterator.hasNext()) {
			numOfLatticePoints++;			
			
			Xn.assignTimes( reB, n );
			Yn.assignTimes( imB, n );
			
			final double nXn = n.dot(Xn);
			final double nYn = n.dot(Yn);
			
			final double nx = n.dot(x);
			final double ny = n.dot(y);
			
			term.assignExp( nXn/2 + nx, nYn/2 +ny );
			
			ComplexVector.dotBilinear( n, X, nX);
			ComplexVector.dotBilinear( n, Y, nY);
			
			termX.assignTimes( term, nX);
			termY.assignTimes( term, nY);
			
			termXY.assignTimes( termX, nY );
			
			thetaSumZ.assignPlus( term );
			thetaSumX.assignPlus( termX );
			thetaSumY.assignPlus( termY );
			
			thetaSumXY.assignPlus( termXY );
		}
	}
	
	final void ddThetaSum(
		final ComplexVector Z,
		final ComplexVector X,
		final ComplexVector Y,
		final Complex thetaSumZ,
		final Complex thetaSumX,
		final Complex thetaSumY,
		final Complex thetaSumXY) {
		
		if( this.uniformApproximation ) {
			ddThetaSumUniform( Z, X, Y, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY );
		} else {
			ddThetaSumPointwise( Z, X, Y, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY );
		}
	}
	
	final void ddThetaSumUniform(
			final ComplexVector Z,
			final ComplexVector X,
			final ComplexVector Y,
			final Complex thetaSumZ,
			final Complex thetaSumX,
			final Complex thetaSumY,
			final Complex thetaSumXY) {
		
		final double[] zRe = Z.re, zIm = Z.im;
		final double[] uRe = X.re, uIm = X.im;
		final double[] vRe = Y.re, vIm = Y.im;

		final double[][] intLatticePointsRe = latticePoints.re;

		final double[] expOfHalfBnnRe = expOfHalfBnn.re, expOfHalfBnnIm = expOfHalfBnn.im;

		thetaSumZ.assign(1, 0);
		thetaSumX.assign(0, 0);
		thetaSumY.assign(0, 0);
		thetaSumXY.assign(0, 0);

		for (int i = 1; i < numOfLatticePoints; i++) {

			final double[] nRe = intLatticePointsRe[i];

			double nZRe = 0, nZIm = 0;
			double nXRe = 0, nXIm = 0;
			double nYRe = 0, nYIm = 0;

			for (int j = 0; j < dim; j++) {
				final double n = nRe[j];

				nZRe += zRe[j] * n;
				nZIm += zIm[j] * n;

				nXRe += uRe[j] * n;
				nXIm += uIm[j] * n;

				nYRe += vRe[j] * n;
				nYIm += vIm[j] * n;
			}

			expOfNZ.assignExp(nZRe, nZIm);
			invOfExpOfNZ.assignInvert(expOfNZ);

			// compute thetaSumZ
			term.assignPlus(expOfNZ, invOfExpOfNZ);
			term.assignTimes(expOfHalfBnn.re[i], expOfHalfBnn.im[i]);
			thetaSumZ.assignPlus(term);

			// compute thetaSumXY
			term.assignTimes(nXRe, nXIm);
			term.assignTimes(nYRe, nYIm);
			thetaSumXY.assignPlus(term);

			// compute thetaSumX and thetaSumY
			term.assignMinus(expOfNZ, invOfExpOfNZ);
			term.assignTimes(expOfHalfBnn.re[i], expOfHalfBnn.im[i]);

			tmp.assign(term);

			term.assignTimes(nXRe, nXIm);
			thetaSumX.assignPlus(term);

			term.assign(tmp);

			term.assignTimes(nYRe, nYIm);
			thetaSumY.assignPlus(term);
		}
	}

	private final void getCombinedFactor(Complex f) {

		f.assignPlus(modular.factor, transform.factor);
	}

	private final void getDerivativeOfCombinedFactor(ComplexVector X, Complex dFByDX) {

		tmpVector.assignTimes(transform.M, modular.H);

		tmpVector.assignPlus(modular.R);
		tmpVector.assignMinus(modular.ATrZ);
		tmpVector.assignMinus(modular.AZ);

		ComplexVector.dotBilinear(tmpVector, X, dFByDX);
	}

	private final void getDerivativeOfCombinedFactor(ComplexVector X, ComplexVector Y, Complex ddFByDXDY) {

		tmpVector.assignTimes(modular.A, X);

		ComplexVector.dotBilinear(tmpVector, Y, ddFByDXDY);

		tmpVector.assignTimes(modular.A, Y);

		ComplexVector.dotBilinear(tmpVector, X, tmp);

		ddFByDXDY.assignPlus(tmp);
		ddFByDXDY.assignNeg();
	}

	private final void getDervativeOfCombinedZTransformation(ComplexVector X, ComplexVector dTByDX) {

		dTByDX.assignTimes(modular.H, X);
	}

	private final void getContinousFactor(final Complex factor, final Complex continousCorrection) {

		double continousFactor = -0.5 * transform.zRe.dot(transform.reBInvReZ);

		continousCorrection.assign(transform.factor);
		continousCorrection.re -= continousFactor;

		factor.assignMinus(continousCorrection);

		continousCorrection.assignExp();
	}

	/**
	 * Evaluates the Riemann theta function at <code>Z</code>.
	 * <p align=center>
	 *   <code>
	 *     &theta;(z|B) = exp( factor ) &middot; thetaSumZ
	 *   </code>
	 * </p>
	 * @param Z argument vector
	 * @param factor exponential part
	 * @param thetaSumZ oscillatory part
	 */
	public final void theta(final ComplexVector Z, final Complex factor, final Complex thetaSumZ) {

		if (modularIsId) {

			transform.setZ(Z);
			factor.assign(transform.factor);

		} else {

			modular.setZ(Z);
			transform.setZ(modular.tZ);
			getCombinedFactor(factor);
		}

		thetaSum(transform.transfromedZ, thetaSumZ);

		getContinousFactor(factor, continousCorrection);

		thetaSumZ.assignTimes(continousCorrection);
	}

	/**
	 * Evaluates the Riemann theta function and its first derivative
	 * in the <code>X</code> direction at <code>Z</code>.
	 * <p align=center>
	 * <code>
	 * <table>
	 * <tr> <td>              &theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumZ</td> </tr>
	 * <tr> <td> D<sub>X</sub>&theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumX</td> </tr>
	 * </table>
	 * </code>
	 * <p>
	 * Evaluating the Riemann theta funciton and its first derivative simultanesly
	 * can be performed with almost no extra cost.
	 * @param Z argument vector
	 * @param X direction of derivative
	 * @param factor exponential part
	 * @param thetaSumZ oscillatory part
	 * @param thetaSumX oscillatory part of first derivative in <code>X</code> direction
	 */
	public final void dTheta(final ComplexVector Z, final ComplexVector X, final Complex factor, final Complex thetaSumZ, final Complex thetaSumX) {

		if (modularIsId) {

			transform.setZ(Z);

			factor.assign(transform.factor);

			dThetaSum(transform.transfromedZ, X, thetaSumZ, thetaSumX);

			ComplexVector.dotBilinear(X, transform.M, MX);

			tmp.assignTimes(MX, thetaSumZ);
			thetaSumX.assignPlus(tmp);

		} else {

			modular.setZ(Z);
			transform.setZ(modular.tZ);

			getCombinedFactor(factor);

			getDerivativeOfCombinedFactor(X, dFactorByDX);
			getDervativeOfCombinedZTransformation(X, dTByDX);

			dThetaSum(transform.transfromedZ, dTByDX, thetaSumZ, thetaSumX);

			tmp.assignTimes(dFactorByDX, thetaSumZ);
			thetaSumX.assignPlus(tmp);
		}

		getContinousFactor(factor, continousCorrection);

		thetaSumZ.assignTimes(continousCorrection);
		thetaSumX.assignTimes(continousCorrection);
	}

	/**
	 * Evaluates the Riemann theta function, its first derivatives
	 * in the <code>X</code> and <code>Y</code> direction and,
	 * its second derivative into the same direction at <code>Z</code>.
	 * <p align=center>
	 * <code>
	 * <table>
	 * <tr> <td>                       &theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumZ </td> </tr>
	 * <tr> <td>          D<sub>X</sub>&theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumX </td> </tr>
	 * <tr> <td>          D<sub>Y</sub>&theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumY </td> </tr>
	 * <tr> <td>  D&sup2;<sub>X,Y</sub>&theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumXY</td> </tr>
	 * </table>
	 * </code>
	 * <p>
	 * Evaluating the Riemann theta funciton and its first derivative simultanesly
	 * can be performed with almost no extra cost.
	 * @param Z argument vector
	 * @param X direction of derivative
	 * @param Y direction of derivative
	 * @param factor exponential part
	 * @param thetaSumZ oscillatory part
	 * @param thetaSumX oscillatory part of first derivative in <code>X</code> direction
	 * @param thetaSumY oscillatory part of first derivative in <code>Y</code> direction
	 * @param thetaSumXY oscillatory part of second derivative
	 */
	public final void ddTheta(
		final ComplexVector Z,
		final ComplexVector X,
		final ComplexVector Y,
		final Complex factor,
		final Complex thetaSumZ,
		final Complex thetaSumX,
		final Complex thetaSumY,
		final Complex thetaSumXY) {

		if (modularIsId) {

			transform.setZ(Z);

			factor.assign(transform.factor);

			ddThetaSum(transform.transfromedZ, X, Y, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY);

			ComplexVector.dotBilinear(X, transform.M, MX);
			ComplexVector.dotBilinear(Y, transform.M, MY);

			// compute thetaSumXY
			tmp.assignTimes(MX, thetaSumY);
			thetaSumXY.assignPlus(tmp);
			tmp.assignTimes(MY, thetaSumX);
			thetaSumXY.assignPlus(tmp);
			tmp.assignTimes(MY, thetaSumZ);
			tmp.assignTimes(MX);
			thetaSumXY.assignPlus(tmp);

			// compute thetaSumX
			tmp.assignTimes(MX, thetaSumZ);
			thetaSumX.assignPlus(tmp);
			// compute thetaSumY
			tmp.assignTimes(MY, thetaSumZ);
			thetaSumY.assignPlus(tmp);

		} else {

			modular.setZ(Z);
			transform.setZ(modular.tZ);

			getCombinedFactor(factor);

			getDerivativeOfCombinedFactor(X, dFactorByDX);
			getDerivativeOfCombinedFactor(Y, dFactorByDY);

			getDerivativeOfCombinedFactor(X, Y, ddFactorByDXDY);

			getDervativeOfCombinedZTransformation(X, dTByDX);
			getDervativeOfCombinedZTransformation(Y, dTByDY);

			ddThetaSum(transform.transfromedZ, dTByDX, dTByDY, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY);

			// compute thetaSumXY
			tmp.assignTimes(dFactorByDX, dFactorByDY);
			tmp.assignPlus(ddFactorByDXDY);
			tmp.assignTimes(thetaSumZ);
			thetaSumXY.assignPlus(tmp);

			tmp.assignTimes(dFactorByDY, thetaSumX);
			thetaSumXY.assignPlus(tmp);

			tmp.assignTimes(dFactorByDX, thetaSumY);
			thetaSumXY.assignPlus(tmp);

			// compute thetaSumX
			tmp.assignTimes(dFactorByDX, thetaSumZ);
			thetaSumX.assignPlus(tmp);

			// compute thetaSumY
			tmp.assignTimes(dFactorByDY, thetaSumZ);
			thetaSumY.assignPlus(tmp);

		}

		getContinousFactor(factor, continousCorrection);

		thetaSumZ.assignTimes(continousCorrection);
		thetaSumX.assignTimes(continousCorrection);
		thetaSumY.assignTimes(continousCorrection);
		thetaSumXY.assignTimes(continousCorrection);
	}

}
