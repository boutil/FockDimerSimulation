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
import de.jtem.mfc.field.Complex;

/**
 * Riemann theta function with characteristics.
 * <p align=center>
 *   <code>
 *     &theta;[&alpha;,&beta;](z|B) =  exp( &frac12; ( z + &pi;i&beta; + &frac14;B&alpha; , &alpha; ) )
 *                                     &theta;( z + &pi;i&beta; + &frac12;B&alpha |B)
 *   </code>,
 * </p>
 * where
 * <code>
 * z<small>&isin;</small>C<sup>g</sup>,
 * &alpha;,&beta;<small>&isin;</small>Z<sup>g</sup>
 * </code> and
 * <code>
 * B<small>&isin;</small>C<sup>g<small>x</small>g</sup>
 </code>
 * where <code>B</code> is symmetric and has strictly negative definite real part.
 * <p>
 * The pair <code>[&alpha;,&beta;]</code> is called characteristic. As follows
 * from the definition of the Riemann theta function with characteristics,
 * that there are only <code>4<sup>dim</sup></code> different characteristics which result
 * in different Riemann theta functions. These can be divided into odd and even ones. A
 * Characteristic is called even, iff
 * <p align=center>
 *   <code>
 *     ( &alpha;,&beta; ) = 0 mod 2
 *   </code>,
 * </p>
 * and odd otherwise. There are <code>2<sup>dim-1</sup>( 2<sup>dim-1</sup> + 1 )</code> even and
 * <code>2<sup>dim-1</sup>( 2<sup>dim-1</sup> - 1 )</code> odd characteristics.
 * <p>
 * Each instance of <code>ThetaWithChar</code> is associated to an
 * instance of <code>Theta</code>. These associated Riemann theta function
 * can be shared by different Riemann theta functions with characteristics.
 * Suppose you want to evaluate
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td> &theta;[&alpha,&beta](z+a|B)    </td> </tr>
 * <tr> <td> <hr size=1 noshade>             </td> </tr>
 * <tr> <td> &theta;[&gamma,&delta](z+b|B) . </td> </tr>
 * </table>
 * </code>
 * <p>
 * The most efficient way to implement this is to create two instances of
 * <code>ThetaWithChar</code> sharing one instance of <code>Theta</code>:
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td> Theta associatedTheta = new Theta( B, tol );                         </td></tr>
 * <tr> <td>                                                                      </td></tr>
 * <tr> <td> ThetaWithChar thetaAlphaBeta  = new ThetaWithChar( associatedTheta );</td></tr>
 * <tr> <td> ThetaWithChar thetaGammaDelta = new ThetaWithChar( associatedTheta );</td></tr>
 * <tr> <td>                                                                      </td></tr>
 * <tr> <td> thetaAlphaBeta .setAlpha( alpha );                                   </td></tr>
 * <tr> <td> thetaAlphaBeta .setBeta(  beta  );                                   </td></tr>
 * <tr> <td>                                                                      </td></tr>
 * <tr> <td> thetaGammaDelta.setAlpha( gamma );                                   </td></tr>
 * <tr> <td> thetaGammaDelta.setBeta(  delta );                                   </td></tr>
 * </table>
 * </code>
 * <p>
 * The benefit of this is that you only have to initialze one Riemann theta function, and
 * the initialization can be several orders more &quot;expensive&quot; than its evaluation.
 * But be aware that the three objects <code>thetaAlphaBeta, thetaGammaDelta</code>, and
 * <code>associatedTheta</code> are dependend.
 *
 * @see Theta


 * @author Markus Schmies
 **/
public class ThetaWithChar extends AbstractTheta implements Serializable, Cloneable {

    /**
     * Name conventions in this class:
     *
     *  &theta;[alpha,beta](z,B)
     *     = exp( &frac12; <z,alpha> ) exp( phi ) &theta;( z + tau, B )
     *
     * with   phi = &frac12; < Pi i beta + &frac14; B alpha,alpha>
     * and    tau =       Pi i beta + &frac12; B alpha
     *
     */

    private static final long serialVersionUID = 1L;

    private final Complex tmp              = new Complex();
    private final Complex i2Pi             = new Complex( 0, 2*Math.PI);

    private final ComplexVector tau        = new ComplexVector();
    private final ComplexVector halfAlpha  = new ComplexVector();
    private final ComplexVector halfBeta   = new ComplexVector();
    private final ComplexVector halfBAlpha = new ComplexVector();
    private final ComplexVector tauPlusZ   = new ComplexVector();

    private final Complex phi              = new Complex();

    private final Complex halfXAlpha       = new Complex();
    private final Complex uAlphaHalf       = new Complex();
    private final Complex zAlphaHalf       = new Complex();


    IntegerVector alpha = new IntegerVector ();
    IntegerVector beta  = new IntegerVector ();

    Theta theta;

    int dim;

    int parityOfSpinStructure;
    int numOfCharacteristics;

    long lastUpdate = Long.MIN_VALUE;

    boolean uptodate = false;

    /**
     * Creates a Riemann theta function with characteristics zero
     * and with prescribed <code>periodMatrix</code>.
     * The error tolerance is by default <code>e-7</code>, the fill factor error
     * and Siegel`s reduction is used.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public ThetaWithChar( final ComplexMatrix periodMatrix ) {

	theta = new Theta( periodMatrix );

        update();
    }

    /**
     * Creates a Riemann theta function with characteristics zero
     * and with prescribed <code>periodMatrix</code>
     * and error tolerance <code>tol</code>.
     * Both, fill factor error and Siegel`s reduction is used.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param tol positive number
     */
    public ThetaWithChar( final ComplexMatrix periodMatrix, final double tol ) {

	theta = new Theta( periodMatrix, tol );

        update();
    }

    /**
     * Creates a Riemann theta function with characteristics zero
     * and with prescribed <code>periodMatrix</code>
     * and error tolerance <code>tol</code>.
     * It further enables the user to switch off the Siegel`s Reduction.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param tol positive number
     * @param performSiegelReduction  controles whether or not Siegel`s reduction is performed
     */
    public ThetaWithChar( final ComplexMatrix periodMatrix, final double tol,
			  final boolean performSiegelReduction ) {

	theta = new Theta( periodMatrix, tol, performSiegelReduction );

        update();
    }

    /**
     * Creates a Riemann theta function with characteristics zero and
     * with prescribed <code>periodMatrix</code>
     * and error tolerance <code>tol</code>.
     * It further enables the user to switch off the Siegel`s Reduction and the
     * heuristic of the fill factor error.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param tol positive number
     * @param performSiegelReduction  controles whether or not Siegel`s reduction is performed
     * @param useFillFactorError controles whether the fill factor or the 100% error is used.
     */
    public ThetaWithChar( final ComplexMatrix periodMatrix, final double tol,
			  final boolean performSiegelReduction, final boolean useFillFactorError ) {

	theta = new Theta( periodMatrix, tol, performSiegelReduction, useFillFactorError );

        update();
    }

    /**
     * Creates a Riemann theta function with characteristics zero and
     * with prescribed <code>periodMatrix</code>
     * and error tolerance <code>tol</code>.
     * It further enables the user to switch off the Siegel`s Reduction and the
     * heuristic of the fill factor error.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param tol positive number
     * @param performSiegelReduction  controles whether or not Siegel`s reduction is performed
     * @param useFillFactorError controles whether the fill factor or the 100% error is used.
     * @param uniformApproximation controles whether to use uniform or pointwise approximation.
     */
    public ThetaWithChar( final ComplexMatrix periodMatrix, final double tol,
    		final boolean performSiegelReduction, final boolean useFillFactorError, final boolean uniformApproximation ) {

    	theta = new Theta( periodMatrix, tol, performSiegelReduction, useFillFactorError, uniformApproximation );

    	update();
    }
    
    /**
     * Creates a Riemann theta function with characteristics zero which
     * is associated to <code>theta</code>.
     * Associated means that <code>this</code> will deligate
     * most of its functionallity to <code>theta</code>. Thus all other Riemann
     * theta funciton objects with characteristics
     * which are also associated to <code>theta</code>
     * are implicitly associated with <code>this</code>.
     * @param theta associated Riemann theta function
     */
    public ThetaWithChar( final Theta theta ) {
	this.theta = theta;
	uptodate = false;

        update();
    }

    /**
     * Creates a Riemann theta function with characteristics zero which
     * shares with <code>thetaWithChar</code> its associated Riemann theta function.
     * @param thetaWithChar associated Riemann theta function with characteristics
     */
    public ThetaWithChar( final ThetaWithChar thetaWithChar ) {
	this( thetaWithChar.theta );
    }

    /**
     * Returns the associated Riemann theta function.
     */
    public final Theta getAssociatedTheta() {
	return theta;
    }
    /**
     * Returns the fill factor.
     */
    public final double getFillFactor() {
	return theta.getFillFactor();
    }

    /**
     * Returns whether the heuristic of the fill factor error is used or not.
     */
    public final boolean getUseFillFactorError() {
	return theta.getUseFillFactorError();
    }

    /**
     * Switches the heuristic fill factor error on or off.
     * Affects the associated Riemann theta function and therefor
     * all other instances of <code>ThetaWithChar</code> which
     * share the same associated Riemann theta function.
     * @param useFillFactorError controles whether the fill factor or the 100% error is used.
     */
    public final void setUseFillFactorError( final boolean useFillFactorError ) {
	theta.setUseFillFactorError( useFillFactorError );
    }

    /**
     * @deprecated
     * Returns the error tolerance for the oscillatory part of the theta function
     * with characteristics.
     */
    public final double getTol() {
	return theta.getTol();
    }

    /**
     * @deprecated
     * Sets the error tolerance for the oscillatory part of the theta function.
     * Affects the associated Riemann theta function and therefor
     * all other instances of <code>ThetaWithChar</code> which
     * share the same associated Riemann theta function.
     * @param tol positive number
     */
    public final void setTol( final double tol ) {
	theta.setTol( tol );
    }

    /**
     * Returns period matrix on which the Riemann theta function with characteristics operates.
     */
    public final ComplexMatrix getPeriodMatrix() {
	return theta.getPeriodMatrix();
    }

    /**
     * Sets periodMatrix on which the Riemann theta function with characteristics operates.
     * Affects the associated Riemann theta function and therefor
     * all other instances of <code>ThetaWithChar</code> which
     * share the same associated Riemann theta function.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public final void setPeriodMatrix( final ComplexMatrix periodMatrix ) {

	theta.setPeriodMatrix( periodMatrix );
    }

    /**
     * Returns number of characteristics.
     * @return <code>4<sup>dim</sup></code>
     */
    public final int getNumOfCharacteristics() {
	update();
	return numOfCharacteristics;
    }

    private void updateNumOfCharacterstics() {
	numOfCharacteristics = 1;
	for( int i=0; i<dim; i++)
	    numOfCharacteristics  *= 4;
    }

    /**
     * Returns number of even characteristics.
     * @return <code>2<sup>dim-1</sup>( 2<sup>dim-1</sup> + 1 )</code>
     */
    public final int getNumOfEvenCharacteristics() {
	update();
	int num = 1;
	for( int i=1; i<dim; i++)
	    num  *= 2;
	return num * ( 2 * num + 1 );
    }

    /**
     * Returns number of odd characteristics.
     * @return <code>2<sup>dim-1</sup>( 2<sup>dim-1</sup> - 1 )</code>
     */
    public final int getNumOfOddCharacteristics() {
	update();
	int num = 1;
	for( int i=1; i<dim; i++)
	    num  *= 2;
	return num * ( 2 * num - 1 );
    }

    final void update() {

	if( uptodate && lastUpdate > theta.lastChangeOfPeriodMatrix )
	    return;

	setDim( theta.getDim() );

	parityOfSpinStructure = alpha.dot( beta );

	for( int i=0; i<dim; i++ ) {
	  halfAlpha.re[i] = alpha.re[i] / 2.;
	  halfBeta. re[i] = beta. re[i] / 2.;
	}

	halfBAlpha.assignTimes( theta.periodMatrix, halfAlpha );

	tau.assignTimes( halfBeta, i2Pi );
	tau.assignPlus( halfBAlpha );

	ComplexVector.dotBilinear( halfBAlpha, halfAlpha, phi );
	ComplexVector.dotBilinear( halfBeta,   halfAlpha, tmp );

	tmp.assignTimes( i2Pi );
	phi.assignTimes( 0.5 );
	phi.assignPlus( tmp );

	uptodate = true;

	lastUpdate = System.currentTimeMillis();
    }

    /**
     * Returns the <code>i<sup>th</sup></code>-component of
     * the <code>alpha</code> charateristic.
     * @param ith non negative, smaller than <code>dim</code>
     * @return Returns the <code>i<sup>th</sup></code>-component of
     * the <code>alpha</code> charateristic
     */
     public final int getAlpha( final int ith ) {
	 return alpha.get( ith );
     }

    /**
     * Sets the <code>i<sup>th</sup></code>-component of
     * the <code>alpha</code> charateristic with <code>value</code>.
     * @param ith non negative, smaller than <code>dim</code>
     * @param value new value for the <code>i<sup>th</sup></code>-component of
     * the <code>alpha</code> charateristic
     */
    public final void setAlpha( final int ith, final int value ) {
	if( value == getAlpha(ith))
		return;
	alpha = new IntegerVector ( alpha );
	alpha.set( ith, value );
	uptodate = false;
    }

     /**
     * Returns the <code>i<sup>th</sup></code>-component of
     * the <code>beta</code> charateristic.
     * @param ith non negative, smaller than <code>dim</code>
     * @return Returns the <code>i<sup>th</sup></code>-component of
     * the <code>beta</code> charateristic
     */
    public final int getBeta( final int ith ) {
	return beta.get( ith );
    }

    /**
     * Sets the <code>i<sup>th</sup></code>-component of
     * the <code>beta</code> charateristic with <code>value</code>.
     * @param ith non negative, smaller than <code>dim</code>
     * @param value new value for the <code>i<sup>th</sup></code>-component of
     * the <code>beta</code> charateristic
     */
    public final void setBeta( final int ith, final int value ) {
	if( value == getBeta(ith))
			return;
	beta = new IntegerVector ( beta );
	beta.set( ith, value );
	uptodate = false;
    }

    /**
     * Returns <code>alpha</code> charateristic.
     */
    public final IntegerVector getAlpha() {
	return new IntegerVector ( alpha );
    }

    /**
     * Sets the <code>alpha</code> charateristic to <code>alpha</code> .
     * @param alpha integer vector of size <code>dim</code>
     */
    public final void setAlpha( final IntegerVector alpha ) {
	if( alpha.size() != dim )
	    throw new IllegalArgumentException
                ( "charctersics has wrong size; right size is " + dim );
	if( alpha.equals(this.alpha))
		return;
	this.alpha = new IntegerVector ( alpha );
	uptodate = false;
    }

    /**
     * Returns <code>beta</code> charateristic.
     */
    public final IntegerVector getBeta() {
	return new IntegerVector ( beta );
    }

    /**
     * Sets the <code>beta</code> charateristic to <code>beta</code> .
     * @param beta integer vector of size <code>dim</code>
     */
    public final void setBeta( final IntegerVector beta ) {
	if( beta.size() != dim )
	    throw new IllegalArgumentException( "charctersics has wrong size" );

	if( beta.equals(this.beta))
		return;
	this.beta = new IntegerVector ( beta );
	uptodate = false;
    }

    /**
     * @depricated Returns parity of spin structure.
     * @return <alpha,beta>
     */
    public final int getParityOfSpinStructure() {
	update();
	return parityOfSpinStructure;
    }

    /**
     * Retunrs <code>true</code> if <code>this</code> has an even characteristic.
     * @return  <code>(alpha,beta) % 2 == 0</code>
     */
    public final boolean isEvenCharacteristic() {
	return alpha.dot( beta ) % 2 == 0;
    }

    /**
     * Retunrs <code>true</code> if <code>this</code> has an odd characteristic.
     * @return  <code>(alpha,beta) % 2 != 0 </code>
     */
    public final boolean isOddCharacteristic() {
	return alpha.dot( beta ) % 2 != 0;
    }

    /** Returns the dimension of the complex vector space on which the Riemann theta
	functions operates.
    */
    public final int getDim() {
	update();
	return dim;
    }

    final void setDim(int v) {

	if( v != dim ) {

	    dim = v;

	    updateNumOfCharacterstics();

	    theta.setDim( v );

	    alpha.resize( dim );
	    beta.resize( dim );

	    halfAlpha.resize( dim );
	    halfBeta. resize( dim );
	}
    }

    /**
     * Evaluates the Riemann theta function with characteristics at <code>Z</code>.
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
	update();

	/*  &theta;[alpha,beta](z,B)
	 *     = exp( &frac12; <z,alpha> ) exp( phi ) &theta;( z + tau, B )
	 */

	tauPlusZ.assignPlus( Z, tau );

	theta.theta( tauPlusZ, factor, thetaSumZ );

	ComplexVector.dotBilinear( Z, halfAlpha, zAlphaHalf );

	factor.assignPlus( zAlphaHalf );
	factor.assignPlus( phi );
    }



    /**
     * Evaluates the Riemann theta function with characteristics and its first derivative
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
	update();

	/* d &theta;[alpha,beta](z) * dz =
	 *	        exp( phi ) * exp( 1/2 <z alpha> ) *
	 *             ( 1/2 < dz, alpha > &theta;( z + tau ) + d &theta;( z + tau ) dz )
	 */

	tauPlusZ.assignPlus( Z, tau );

	theta.dTheta( tauPlusZ, X, factor, thetaSumZ, thetaSumX );

	ComplexVector.dotBilinear( Z, halfAlpha, zAlphaHalf );
	ComplexVector.dotBilinear( X, halfAlpha, halfXAlpha );

	halfXAlpha.assignTimes( thetaSumZ );
	thetaSumX.assignPlus( halfXAlpha );

	factor.assignPlus( zAlphaHalf );
	factor.assignPlus( phi );
    }

    /**
     * Evaluates the Riemann theta function with characteristics, its first derivatives
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
			       final Complex thetaSumY,
			       final Complex thetaSumX,
			       final Complex thetaSumXY ) {
	update();

	/* dd theta[alpha,beta](z)(Y,X) =
	 *   exp( phi ) * exp( 1/2 <z alpha> ) * (  dd theta(z+tau)(Y,X)
	 *                                        + 1/2 <Y,alpha> d theta(z+tau)(X)
	 *                                        + 1/2 <X,alpha> d theta(z+tau)(Y)
	 *                                        + 1/4<Y,alpha><X,alpha> theta(z+tau)(Y) ) */

	tauPlusZ.assignPlus( Z, tau );

	theta.ddTheta( tauPlusZ, X, Y, factor, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY );

	ComplexVector.dotBilinear( Z, halfAlpha, zAlphaHalf );
	ComplexVector.dotBilinear( X, halfAlpha, halfXAlpha );
	ComplexVector.dotBilinear( Y, halfAlpha, uAlphaHalf );

	tmp.assignTimes( uAlphaHalf, thetaSumX );
	thetaSumXY.assignPlus( tmp );

	tmp.assignTimes( halfXAlpha, thetaSumY );
	thetaSumXY.assignPlus( tmp );

	tmp.assignTimes( halfXAlpha, uAlphaHalf );
	tmp.assignTimes( thetaSumZ );
	thetaSumXY.assignPlus( tmp );

	halfXAlpha.assignTimes( thetaSumZ );
	thetaSumX. assignPlus( halfXAlpha );

	uAlphaHalf.assignTimes( thetaSumZ );
	thetaSumY. assignPlus( uAlphaHalf );

	factor.assignPlus( zAlphaHalf );
	factor.assignPlus( phi );
    }

}
















