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
import de.jtem.blas.IntegerMatrix;
import de.jtem.blas.IntegerVector;
import de.jtem.mfc.field.Complex;

/**
 * Supports the use of the modular property of Riemann theta functions.
 *
 * Unlike the transformation property is the modular transformation property
 * a deep result in the theory of Riemann theta functions. For all elements
 * <p align=center>
 * <code>
 * <table>
 * <tr>
 * <td> &sigma; =
 * </td>
 * <td>
 * <table>
 * <tr> <td>/ a b \</td></tr>
 * <tr> <td>\ c d /</td></tr>
 * </table>
 * </td>
 * </tr>
 * </table>
 * </code>
 * <p>
 * of the {@link ModularGroup} of degree <code>2g</code> holds:
 * <p align=center>
 *   <code>
 *     &theta;(z|B) = k exp( f<sub>&sigma;</sub>(z|B) ) &theta;( &sigma;(z|B) | &sigma;(B) )
 *   </code>,
 * </p>
 * with
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td align=right> f<sub>&sigma;</sub>(z|B)</td> <td> = </td>
 *      <td>R<sub>&sigma;</sub>(B)<sup>t</sup>z -
 *          z<sup>t</sup>A<sub>&sigma;</sub>(B)z +
 *          &delta;<sub>&sigma;</sub>(B)                      </td> </tr>
 * <tr> <td align=right> &sigma;(z|B)</td>             <td> = </td>
 *      <td>H<sub>&sigma;</sub>(B)z + S<sub>&sigma;</sub>(B) </td> </tr>
 * <tr> <td align=right> &sigma;(B)</td>               <td> = </td>
 *      <td>2&pi;i(aB+2&pi;ib)(cB+2&pi;id)<sup>-1</sup></td> </tr>
 * </table>
 * </code>
 * <p> and further with
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td align=right>R<sub>&sigma;</sub>(B)</td>       <td> = </td>
 *      <td>&pi;i (cB + 2&pi;id)<sup>-1</sup>Diag(cd<sup>t</sup>) </td> </tr>
 * <tr> <td align=right>A<sub>&sigma;</sub>(B)</td>       <td> = </td>
 *      <td>&frac12; ((cB + 2&pi;id)<sup>-1</sup>)c</td> </tr>
 * <tr> <td align=right>&delta;<sub>&sigma;</sub>(B)</td> <td> = </td>
 *      <td>&frac12; ( &frac14; Diag(cd<sup>t</sup>)<sup>t</sup>&sigma;(B)Diag(cd<sup>t</sup>) -
 *                    g log(2&pi;i) - 2log( (det(cB +2&pi;id) )<sup>&frac12;</sup> )</td> </tr>
 * <tr> <td align=right>H<sub>&sigma;</sub>(B)</td>       <td> = </td>
 *      <td>2&pi;i( (cB + 2&pi;id)<sup>-1</sup>)<sup>t</sup> </td> </tr>
 * <tr> <td align=right>S<sub>&sigma;</sub>(B)</td>       <td> = </td>
 *      <td>&pi;i Diag(ab<sup>t</sup>) + &frac12; &sigma;(B)Diag(cd<sup>t</sup>)</td> </tr>
 * </table>
 * </code>
 * <p>
 * Finally <code>k</code> is a nasty 8<sup>th</sup>-root of unity,
 * which we do not control. This is usually neglectable because one only deals
 * with ratios of theta functions.
 * <p>
 * To use this class
 * you have to provide a period Matrix <code>B</code> ({@link #setPeriodMatrix}),
 * an element of the {@link ModularGroup} ({@link #setModular}) and, the
 * argument vector <code>z</code> ({@link #setZ}).
 * <code>ModularPropertySupport</code> than determines the transformed argument
 * <code>&sigma;(z|B)</code>, the exponential transformation factor
 * <code>f<sub>&sigma;</sub>(z|B)</code> and,
 * the transformed period matrix <code>&sigma;(B)</code> which can be accessed by calling
 * {@link #getTransformedZ()}, {@link #getFactor()} and, {@link #getTransformedPeriodMatrix()}.
 * Further are
 * <code>
 * R<sub>&sigma;</sub>(B),
 * A<sub>&sigma;</sub>(B),
 * &delta;<sub>&sigma;</sub>(B),
 * H<sub>&sigma;</sub>(B)
 * </code> and, <code>
 * S<sub>&sigma;</sub>(B)
 * </code>
 * attributes of this class.
 * @see SiegelReduction
 * @see ModularTransformation
 * @see TransformPropertySupport
 * @see Theta
 * @author Markus Schmies
 */
public class ModularPropertySupport implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    final static Complex PI1I = new Complex(0, Math.PI);
    final static Complex PI2I = new Complex(0, 2 * Math.PI);
    final static Complex PI4I = new Complex(0, 4 * Math.PI);

    int dim;

    ModularTransformation modular;

    IntegerMatrix a;
    IntegerMatrix aTr;
    IntegerMatrix b;
    IntegerMatrix bTr;
    IntegerMatrix c;
    IntegerMatrix cTr;
    IntegerMatrix d;
    IntegerMatrix dTr;

    ComplexMatrix PM, tPM;
    ComplexVector Z, tZ;

    ComplexVector AZ   = new ComplexVector();
    ComplexVector ATrZ = new ComplexVector();
    ComplexVector HZ   = new ComplexVector();

    Complex factor = new Complex();
    Complex delta  = new Complex();
    IntegerMatrix p = new IntegerMatrix ();

    ComplexVector V = new ComplexVector();
    ComplexVector W = new ComplexVector();
    ComplexVector S = new ComplexVector();
    ComplexVector R = new ComplexVector();

    Complex VZ  = new Complex();
    Complex vV  = new Complex();
    Complex tmp = new Complex();

    ComplexMatrix E = new ComplexMatrix();
    ComplexMatrix F = new ComplexMatrix();
    ComplexMatrix G = new ComplexMatrix();
    ComplexMatrix A = new ComplexMatrix();
    ComplexMatrix H = new ComplexMatrix();

    Complex detOfF;
    Complex logOfSqrtOfDetOfF = new Complex();
    Complex logOfK            = new Complex();

    IntegerVector transformedAlpha;
    IntegerVector transformedBeta;

    /**
     * Creates instance with presrcibed period matrix.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public ModularPropertySupport( final ComplexMatrix periodMatrix ) {
	this( periodMatrix.getNumRows() );
	setPeriodMatrix( periodMatrix );
    }

    /**
     * Creates instance with presrcibed period matrix 
     * and modular transformation <code>&sigma;</code>..
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param sigma element of the modular group
     */
    public ModularPropertySupport( final ComplexMatrix periodMatrix,
    			final ModularTransformation sigma ) {
    	this( periodMatrix.getNumRows() );
    	if(sigma.dim != dim )
    		throw new IllegalArgumentException
			( "dimensions do not match");
    	PM.assign(periodMatrix);
    	modular.assign( sigma );
    	compute();
    }
    
    ModularPropertySupport( final int dim) {
        this.dim = dim;

        modular = new ModularTransformation(dim);

        a = modular.a; aTr = a.transpose();
        b = modular.b; bTr = b.transpose();
        c = modular.c; cTr = c.transpose();
        d = modular.d; dTr = d.transpose();

        PM  = new ComplexMatrix(dim);
        tPM = new ComplexMatrix(dim);

        PM. setDiagonal(2 * Math.PI, 0);
        tPM.setDiagonal(2 * Math.PI, 0);

        Z  = new ComplexVector(dim);
        tZ = new ComplexVector(dim);

        transformedAlpha = new IntegerVector (dim);
        transformedBeta  = new IntegerVector (dim);

        logOfK.assignTimes(PI2I.log(), dim / 2.0);
        logOfK.assignNeg();
    }

    /**
     * Returns period matrix which was transformed.
     */
    public ComplexMatrix getPeriodMatrix() {
        return new ComplexMatrix(PM);
    }

    /**
     * Sets the period matrix which will be transformed.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public void setPeriodMatrix(ComplexMatrix aPeriodMatrix) {
        PM.assign(aPeriodMatrix);
        compute();
    }

    /**
     * Returns transformed period matrix.
     * @return <code>&sigma(B) = 2&pi;i(aB+2&pi;ib)(cB+2&pi;id)<sup>-1</sup></code>
     */
    public ComplexMatrix getTransformedPeriodMatrix() {
        return new ComplexMatrix(tPM);
    }

    /**
     * Returns the argument vector which was transform.
     */
    public ComplexVector getZ() {
        return new ComplexVector(Z);
    }

    /**
     * Sets the argument vector which will be transformed.
     * @param z argument vector
     */
    public void setZ(ComplexVector aZ) {
        Z.assign(aZ);
        AZ.  assignTimes(A, Z);
        ATrZ.assignTimes(Z, A);

        computeTransformedZ();
        computeFactor();
    }

    /**
     * Returns the transformed argument vector.
     * @return <code>&sigma;(z|B) = H<sub>&sigma;</sub>(B)z + S<sub>&sigma;</sub>(B)</code>
     */
    public ComplexVector getTransformedZ() {
        return new ComplexVector(tZ);
    }

    /**
     * Returns the exponential transformation factor.
     * @return <code>f<sub>&sigma;</sub>(z|B) =
     *                       R<sub>&sigma;</sub>(B)<sup>t</sup>z -
     *                       z<sup>t</sup>A<sub>&sigma;</sub>(B)z +
     *                       &delta;<sub>&sigma;</sub>(B)</code>
     */
    public Complex getFactor() {
        return new Complex(factor);
    }

    /**
     * Returns the element of the modular group that determines the transformation.
     */
    public ModularTransformation getModularTransformation() {
        return new ModularTransformation(modular);
    }

    /**
     * Sets the element of the modular group that determines the transformation.
     * @param modular element of the modular group
     */
    public void setModularTransformation( final ModularTransformation modular ) {
    	if(modular.dim != dim )
    		throw new IllegalArgumentException
				( "dimensions do not match");
    	
        this.modular.assign(modular);

        aTr.assignTranspose(a);
        bTr.assignTranspose(b);
        cTr.assignTranspose(c);
        dTr.assignTranspose(d);

        compute();
    }

    /**
     * Returns the constant component of the transformation factor.
     * @return <code>&delta =
     *     &frac12; ( &frac14; Diag(cd<sup>t</sup>)<sup>t</sup>&sigma;(B)Diag(cd<sup>t</sup>) -
     *                    g log(2&pi;i) - 2log( (det(cB +2&pi;id) )<sup>&frac12;</sup> )</code>
     */
    public Complex getDelta() {
	return new Complex( delta );
    }

    /**
     * Returns the linear component of the transformation factor.
     * @return <code>R<sub>&sigma;</sub>(B) =
     *      <td>&pi;i (cB + 2&pi;id)<sup>-1</sup>Diag(cd<sup>t</sup>) </code>
     */
    public ComplexVector getR() {
	return new ComplexVector( R );
    }

    /**
     * Returns the quadratic component of the transformation factor.
     * @return <code>A<sub>&sigma;</sub>(B) =
     *      <td>&frac12; ((cB + 2&pi;id)<sup>-1</sup>)c</code>
     */
    public ComplexMatrix getA() {
	return new ComplexMatrix( A );
    }

    /**
     * Returns the constant component of the argument transformation.
     * @return <code>S<sub>&sigma;</sub>(B) =
     *      <td>&pi;i Diag(ab<sup>t</sup>) + &frac12; &sigma;(B)Diag(cd<sup>t</sup>)</code>
     */
    public ComplexVector getS() {
	return new ComplexVector( S );
    }

    /**
     * Returns the linear component of the argument transformation.
     * @return <code>H<sub>&sigma;</sub>(B) =
     *      <td>2&pi;i( (cB + 2&pi;id)<sup>-1</sup>)<sup>t</sup></code>
     */
    public ComplexMatrix getH() {
	return new ComplexMatrix( H );
    }

    void assignId() {
        modular.assignId();

        aTr.assignTranspose(a);
        bTr.assignTranspose(b);
        cTr.assignTranspose(c);
        dTr.assignTranspose(d);

        compute();
    }

    boolean isId() {
	return modular.isId();
    }

    void computeTransformedChar() {
        p.assignTimes(c, dTr);
        p.getDiagonal(transformedAlpha);

        p.assignTimes(a, bTr);
        p.getDiagonal(transformedBeta);

        V.assignTimes(transformedBeta, PI1I);
        W.assignTimes(tPM, transformedAlpha);
        W.assignDivide(2);

        S.assignPlus(V, W);

        R.assignTimes(transformedAlpha, H);
        R.assignDivide(2);

        // compute delta
        ComplexVector.dotBilinear(W, transformedAlpha, delta);
        delta.assignDivide(4);
        delta.assignMinus(logOfK);
        delta.assignMinus(logOfSqrtOfDetOfF);
    }

    void computeTransformedZ() {
        HZ.assignTimes(H, Z);
        tZ.assignPlus(HZ, S);
    }

    void computeFactor() {
        ComplexVector.dotBilinear(R, Z, factor);

        V.assignTimes(A, Z);

        ComplexVector.dotBilinear(V, Z, VZ);

        factor.assignMinus(VZ);
        factor.assignPlus(delta);
    }

    void compute() {
        H.assignTimes(b, PI2I); // using H temporarily
        E.assignTimes(a, PM);
        E.assignPlus(H);

        H.assignTimes(d, PI2I); // using H temporarily
        F.assignTimes(c, PM);
        F.assignPlus(H);

        detOfF = F.determinant();

        logOfSqrtOfDetOfF.assignSqrt(detOfF);
        logOfSqrtOfDetOfF.assignLog();

        G.assignInvert(F);

        A.assignTimes(G, c);
        A.assignDivide(2);

        AZ.assignTimes(A, Z);

        H.assignTimes(G, PI2I);
        H.assignTranspose();

        tPM.assignTimes(E, G);
        tPM.assignTimes(PI2I);

        computeTransformedChar();
        computeTransformedZ();
        computeFactor();
    }
    
    /**
     * Transforms the period matrix <code>B</code> using
     * the element of the modular group <code>&sigma;</code>.
     * The period matrix is complex a
     * symmetric complex matrix with negative definite real part.
     * <code>&sigma;</code> is of the form:
     * <td>
     * <table>
     * <tr> <td>/ A B \</td><td> </td></tr>
     * <tr> <td>\ C D /</td><td>,</td></tr>
     * </table>
     * </td>
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param sigma element of the modular group
     * @return  2&pi;i(aB+2&pi;ib)(cB+2&pi;id)<sup>-1</sup>
     */
    public static ComplexMatrix transformPeriodMatrix
		( ComplexMatrix B, ModularTransformation sigma ) {
    	return new ModularPropertySupport( B, sigma ).getTransformedPeriodMatrix();
    }
}
