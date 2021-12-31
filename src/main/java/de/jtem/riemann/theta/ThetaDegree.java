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
import de.jtem.blas.IntegerVector;
import de.jtem.blas.RealMatrix;
import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.algebra.linear.decompose.Cholesky;
import de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;
import de.jtem.numericalMethods.geometry.latticeReduction.LLL;

/**
 * Riemann theta function with degree
 * <p>
 * The de.jtem.riemann theta function of degree <code>d</code> and subdegree <code>s</code> is defined by its
 * mulit-dimensional Fourier series
 * <p align=center>
 *   <code>
 *     &theta;<sub>s,d</sub>(z|B) = &sum;<sub>n<small>&isin;</small>Z<sup>g</sup></sub>
 *                                  exp( &frac12;(B(n+s/d),dn+s)+(z,dn+s) )
 *   </code>,
 * </p>
 * where
 * <code>z<small>&isin;</small>C<sup>g</sup>,B<small>&isin;</small>C<sup>g<small>x</small>g</sup></code>,
 * and s<sub>i</sub>&lt;d.
 * <code>B</code> is symmetric and has strictly negative definite real part.
 * <p>
 * The exponential growth in <code>z</code> of the Riemann theta function can be factored
 * out and the remaining oscillatory part can be approximated with a prescribed error
 * tolerance. Thus we have:
 * <p align=center>
 *   <code>
 *     &theta;<sub>s,d</sub>(z|B) = exp( factor(z|B) ) &middot; thetaSum<sub>s,d</sub>(z|B)
 *   </code>,
 * </p>
 * A detailed description can be found in Narasimhan`s book <em>Compact Riemann Surfaces</em>.
 * <p>
 * Separating the exponential factor is crucial when dealing with theta functions. In practice,
 * they often appear as ratios, in which case the exponential growth can be canceled.
 * <p>
 * Evaluating the theta function may easily lead to overflows or numerical
 * instabilities, which should be avoid.
 * <p>
 * To create a theta function object you must provide degree, subdegree and period matrix
 * and you may provide an error tolerance.
 * <p>
 * You may switch off the use of a fill factor error. The fill factor error is a heuristic to
 * improve the accuracy of the error prediction. As mentioned earlier the estimates which lead
 * to the approximation error are rather coarse, which is improved by this heuristic for
 * the price of loosing the certainty of the 100% error. The smaller
 * the fill factor, the bigger its influence on the approximation, which
 * anyway is usually irrelevant for small genus (<code>g<4</code>).
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
 * <tr><td>-</td><td>R.Narasimhan.<em>Compact Riemann Surfaces.</em>
 * Basel, 1992. </td></tr>
 * </table>
 * @see Theta
 * @see TransformPropertySupport
 * @see AbelianFunction
 * @author Markus Schmies/Ulrich Heller
 */

public class ThetaDegree extends AbstractTheta implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    private static final Complex PI2I = new Complex( 0, 2 * Math.PI );

    double tol = 1E-7;

    int dim;

    int degree = 1;

    IntegerVector subDegree = new IntegerVector ();

    ComplexMatrix periodMatrix = new ComplexMatrix();

    TransformPropertySupport transform;

    boolean useFillFactorError     = true;

    long lastChangeOfPeriodMatrix = System.currentTimeMillis();

    LatticePointsForUniformApproximation latticePointsforUniformApproximation;

    RealMatrix latticePoints;

    int numOfLatticePoints;

    ComplexMatrix B;

    RealMatrix reB = new RealMatrix ();
    RealMatrix imB = new RealMatrix ();
    RealMatrix T   = new RealMatrix ();

    ComplexVector halfBnn;

    double radius = 1;
    double lSLV;             // length of shortest lattice vector
    double fillFactor;
    double detOfLattice;
    double lSLV2PowOfDim;    // = Math.pow( lSLV, dim )
    double two2PowOfDim;     // = Math.pow( 2,    dim )
    double gammaOfHalfDim;   // = Gamma.of( dim / 2.0 )
    double sqrtOfPI2PowOfDim; // = Math.pow( Math.sqrt( Math.PI ), dim )

    int highestOrder = 0;

    ErrorRadiusSolverFunction ersf = new ErrorRadiusSolverFunction();

    /* the following latticePointss are used as temporary variables in
       varias methods which makes the whole class not thread save!!! */

    private final Complex tmp                 = new Complex();

    private final Complex continousCorrection = new Complex();
    private final Complex exponent            = new Complex();
    private final Complex expOfExponent       = new Complex();
    private final Complex term                = new Complex();
    private final Complex nX                  = new Complex();
    private final Complex nY                  = new Complex();
    private final Complex nZ                  = new Complex();

    private final Complex MX                  = new Complex();
    private final Complex MY                  = new Complex();

    private final Complex factor              = new Complex();
    private final Complex dFactorByDX         = new Complex();
    private final Complex dFactorByDY         = new Complex();
    private final Complex ddFactorByDXDY      = new Complex();
    private final ComplexVector dTByDX        = new ComplexVector();
    private final ComplexVector dTByDY        = new ComplexVector();
    private final ComplexVector ddTByDXDY     = new ComplexVector();

    private final ComplexVector tmpVector     = new ComplexVector();

    /**
     * Creates a Riemann theta function of degree <code>degree</code> and subdegree
     * <code>subDegree</code> with prescribed <code>periodMatrix</code>.
     * The error tolerance is by default <code>e-7</code>, the fill factor error
     * is used.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param degree is a positve integer
     * @param subDegree is an integer vector with the same dimension as <code>periodMatrix</code>
     */
    public ThetaDegree( final ComplexMatrix periodMatrix, final int degree, final IntegerVector subDegree ) {

        this( periodMatrix, 1e-7, degree, subDegree );
    }

    /**
     * Creates a Riemann theta function of degree <code>degree</code> and subdegree
     * <code>subDegree</code> with prescribed <code>periodMatrix</code> and
     * error tolerance <code>tol</code>, the fill factor error
     * is used.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param degree is a positve integer
     * @param subDegree is an integer vector with the same dimension as <code>periodMatrix</code>
     */
    public ThetaDegree( final ComplexMatrix periodMatrix,
			   final double tol, final int degree, final IntegerVector subDegree ) {
	this.degree=degree;
	this.tol = tol;

	if( subDegree.size() != periodMatrix.getNumCols() )
	  throw new IllegalArgumentException("subDegree has wrong dimension");

	this.subDegree.assign( subDegree );

        setPeriodMatrix( periodMatrix );
    }


    /**
     * Sets the subdegree.
     * This is an integer vector with the same dimension as the period matrix. Its entries
     * must be nonnegative integers smaller than the given degree.
     */
    public void setSubDegree( IntegerVector subDegree ) {

	if( subDegree.size() != periodMatrix.getNumCols() )
	  throw new IllegalArgumentException("subDegree has wrong dimension");

	this.subDegree.assign( subDegree );
	computeLatticePoints();
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
    public final void setUseFillFactorError( final boolean useFillFactorError ) {
	if( useFillFactorError == this.useFillFactorError )
	    return;

	this.useFillFactorError = useFillFactorError;

	computeRadius();
	computeLatticePoints();
    }

    /**
     * Returns the error tolerance for the oscillatory part of the theta function.
     */
    public final double getTol() {
	return tol;
    }

    /**
     * Sets the error tolerance for the oscillatory part of the theta function.
     * @param tol positive number
     */
    public final void setTol( final double tol ) {
	if( this.tol == tol )
	    return;

        this.tol = tol;

	computeRadius();
	computeLatticePoints();
    }

    /** Returns the dimension of the complex vector space on which the Riemann theta
	functions operates.
    */
    public final int getDim() {
	return dim;
    }

    void setDim( final int v ) {

	if( v != dim ) {

	    dim = v;

	    transform = new TransformPropertySupport( dim );

	    gammaOfHalfDim    = Gamma.gamma( dim / 2.0 );
	    two2PowOfDim      = Math.pow( 2, dim );
	    sqrtOfPI2PowOfDim = Math.pow( Math.sqrt( Math.PI ), dim );
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
	return new ComplexMatrix( periodMatrix );
    }

    /**
     * Sets periodMatrix on which the Riemann theta function operates.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public final void setPeriodMatrix( final ComplexMatrix periodMatrix )
    {
	if( this.periodMatrix != null && periodMatrix.equals( this.periodMatrix ) )
	    return;
	if( !periodMatrix.isSquared() )
	    throw new IllegalArgumentException("matrix is not squared");
	if( !periodMatrix.isSymmetric() )
	    throw new IllegalArgumentException("matrix is not symmetric");

	setDim( periodMatrix.getNumRows() );

	this.periodMatrix.assign( periodMatrix );

	lastChangeOfPeriodMatrix = System.currentTimeMillis();

	compute();
    }

    final void compute() {

      B = periodMatrix;

      B.getRe( T );

      T.assignTimes( -1 );

      Cholesky.decompose( T.re, T.re );

      // the lattice:  sqrt(1/2) T * Z^g
      detOfLattice = T.get(0,0);
      for( int i=1; i<dim; i++ )
	detOfLattice *= T.get(i,i);
      detOfLattice /= Math.sqrt( two2PowOfDim );

      T.assignTranspose();

      LLL.reduce( T.re, (int [][])null );

      lSLV = 0;
      for( int i=0; i<dim; i++ )
	lSLV += T.get(0,i) * T.get(0,i);
      lSLV = Math.sqrt( lSLV /= 2);


      lSLV2PowOfDim = Math.pow( lSLV, dim );

      fillFactor = 2 * sqrtOfPI2PowOfDim * lSLV2PowOfDim
	/ two2PowOfDim / dim / gammaOfHalfDim / detOfLattice;

      B.getRe( reB );
      B.getIm( imB );

      B.print( "B" );
      transform.setPeriodMatrix( B );

      // the lattice:  sqrt(1/2) T * Z^g
      if( latticePointsforUniformApproximation == null ) {
	latticePointsforUniformApproximation
	  = new LatticePointsForUniformApproximation( reB.times( -0.5 ) );
      } else{
	latticePointsforUniformApproximation.B.assignTimes( reB, -0.5 );
	latticePointsforUniformApproximation.uptodate = false;
      }

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

System.out.print(degree);
subDegree.print(" -> hier");

	latticePointsforUniformApproximation.setRadius( radius );

	latticePointsforUniformApproximation.update();

	latticePoints = latticePointsforUniformApproximation.latticePoints;
	numOfLatticePoints = latticePointsforUniformApproximation.numOfLatticePoints;

	if( halfBnn == null || halfBnn.size() < numOfLatticePoints) {
	    halfBnn = new ComplexVector( numOfLatticePoints);
	}

	final double [][] BRe = B.re;
	final double [][] BIm = B.im;

	for(int i = 0 ; i < numOfLatticePoints; i++) {

	    final double [] nRe = latticePoints.re[i];

	    double re = 0;
	    double im = 0;


	    for( int j = 0; j < dim; j++ ) {

		final double [] rowBRe = BRe[j];
		final double [] rowBIm = BIm[j];

	        final double leftShift = nRe[j]*degree + subDegree.re[j];
		re += rowBRe[j] * leftShift * leftShift / (2.0 * degree);
		im += rowBIm[j] * leftShift * leftShift / (2.0 * degree);

		for( int k = j+1; k < dim; k++ ) {

	            final double rightShift = nRe[k] + subDegree.re[k]/(double)degree;
		    re += rowBRe[k] * leftShift * rightShift;
		    im += rowBIm[k] * leftShift * rightShift;
		}
	    }

	    halfBnn.set( i, re, im );
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


    final void computeRadius() {

	double harmonicThreshRadius
	    = ( Math.sqrt( dim + 2*highestOrder +
			   Math.sqrt( dim*dim + 8*highestOrder ) ) + lSLV ) / 2;

	try {

	    radius = Math.max( ersf.getRadius( tol, 0 ), harmonicThreshRadius );

	} catch( RuntimeException e ) {

	    radius = harmonicThreshRadius;
	}
    }

    void thetaSum( final ComplexVector Z, final Complex thetaSumZ ) {

	final double [] zRe = Z.re, zIm = Z.im;

	final double [][] intLatticePointsRe = latticePoints.re;

	final double [] halfBnnRe = halfBnn.re;
	final double [] halfBnnIm = halfBnn.im;

	thetaSumZ.assign( 0 );

	for(int i = 0 ; i < numOfLatticePoints; i++){

	    final double [] nRe = intLatticePointsRe[i];

	    double nZRe = 0, nZIm = 0;

	    for( int j = 0; j < dim; j++) {

	        final double shift = nRe[j]*degree + subDegree.re[j];
		nZRe += zRe[j] * shift;
		nZIm += zIm[j] * shift;
	    }

	    exponent.assign( nZRe + halfBnnRe[i], nZIm + halfBnnIm[i] );

	    expOfExponent.assignExp( exponent );

	    thetaSumZ.assignPlus( expOfExponent );
	}
    }

    final void dThetaSum( final ComplexVector Z,
			  final ComplexVector X,
			  final Complex thetaSumZ,
			  final Complex thetaSumX ) {

      throw new RuntimeException( "not implemented");
    }

    final void ddThetaSum( final ComplexVector Z,
			   final ComplexVector X,
			   final ComplexVector Y,
			   final Complex thetaSumZ,
			   final Complex thetaSumX,
			   final Complex thetaSumY,
			   final Complex thetaSumXY ) {
      throw new RuntimeException( "not implemented");
    }

    private final void getContinousFactor( final Complex factor, final Complex continousCorrection ) {

      double continousFactor = -0.5 * transform.zRe.dot( transform.reBInvReZ );

      continousCorrection.assign( transform.factor ); continousCorrection.re -= continousFactor;

      factor.assignMinus( continousCorrection );

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
    public final void theta( final ComplexVector Z, final Complex factor, final Complex thetaSumZ ) {

      transform.setZ( Z );
      factor.assignTimes( transform.factor, degree );

      thetaSum( transform.transfromedZ, thetaSumZ );

      //getContinousFactor( factor, continousCorrection );

      // thetaSumZ.assignTimes( continousCorrection );
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
    public final void dTheta( final ComplexVector Z,
			      final ComplexVector X,
			      final Complex factor,
			      final Complex thetaSumZ,
			      final Complex thetaSumX ) {

      transform.setZ( Z );

      factor.assignTimes( transform.factor, degree );

      dThetaSum( transform.transfromedZ, X, thetaSumZ, thetaSumX );

      ComplexVector.dotBilinear( X, transform.M, MX );

      tmp.assignTimes( MX, thetaSumZ ); thetaSumX.assignPlus( tmp );

      //getContinousFactor( factor, continousCorrection );

      //thetaSumZ.assignTimes( continousCorrection );
      //thetaSumX.assignTimes( continousCorrection );
    }

    /**
     * Evaluates the Riemann theta function, its first derivatives
     * in the <code>X</code> and <code>Y</code> direction and,
     * its second derivative into the same direction at <code>Z</code>.
     * <p align=center>
     * <code>
     * <table>
     * <tr> <td>                      &theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumZ </td> </tr>
     * <tr> <td>         D<sub>X</sub>&theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumX </td> </tr>
     * <tr> <td>         D<sub>Y</sub>&theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumY </td> </tr>
     * <tr> <td> D&sup2;<sub>X,Y</sub>&theta;(z|B) </td> <td> = </td> <td> exp( factor ) &middot; thetaSumXY</td> </tr>
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
    public final void ddTheta( final ComplexVector Z,
			       final ComplexVector X,
			       final ComplexVector Y,
			       final Complex factor,
			       final Complex thetaSumZ,
			       final Complex thetaSumX,
			       final Complex thetaSumY,
			       final Complex thetaSumXY ) {

      transform.setZ( Z );

      factor.assignTimes( transform.factor, degree );

      ddThetaSum( transform.transfromedZ, X, Y, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY );

      ComplexVector.dotBilinear( X, transform.M, MX );
      ComplexVector.dotBilinear( Y, transform.M, MY );

      // compute thetaSumX
      tmp.assignTimes( MX, thetaSumZ ); thetaSumX.assignPlus( tmp );
      // compute thetaSumY
      tmp.assignTimes( MY, thetaSumZ ); thetaSumY.assignPlus( tmp );

      // compute thetaSumXY
      tmp.assignTimes( MX, thetaSumY ); thetaSumXY.assignPlus( tmp );
      tmp.assignTimes( MY, thetaSumX ); thetaSumXY.assignPlus( tmp );
      tmp.assignTimes( MY, thetaSumZ );
      tmp.assignTimes( MX );       thetaSumXY.assignPlus(  tmp );

      //getContinousFactor( factor, continousCorrection );

      //thetaSumZ. assignTimes( continousCorrection );
      //thetaSumX. assignTimes( continousCorrection );
      //thetaSumY. assignTimes( continousCorrection );
      //thetaSumXY.assignTimes( continousCorrection );
    }

}
