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
import de.jtem.blas.RealVector;
import de.jtem.mfc.field.Complex;

/**
 * Supports the use of the transform property of Riemann theta functions.
 *
 * From the definition of the Riemann theta function easily follows, that for
 * all integer vectors <code>M,N&isin;Z<sup>g</sup></code> holds:
 * <p align=center>
 *   <code>
 *     &theta;(z|B) = exp( &frac12; M<sup>t</sup>BM + M<sup>t</sup>z ) &theta;( z + 2&pi;iN + BM | B)
 *   </code>.
 * </p>
 * This property is know as the transform property of Riemann theta functions and it can be used
 * to factor out the exponential growth of the Riemann theta function, by determing the
 * integer vectors <code>M</code> and <code>N</code> such that the absolute value of the
 * transformed argument <code>z+2&pi;iN+BM</code> is as small as possible.
 * This is done by the following vectors:
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td>M(z|B)</td> <td> =  -[</td> <td> re(B)<sup>-1</sup> re(z) ]             </td> </tr>
 * <tr> <td>N(z|B)</td> <td> =  -[</td> <td> &frac12 ( im(z) + im(BM(z|B) )/&pi; ]</td> </tr>
 * </table>
 * </code>
 * <p>
 * To use this class
 * you have to provide a period matrix <code>B</code> ({@link #setPeriodMatrix})
 * and the argument vector <code>z</code> ({@link #setZ}).
 * <code>TransformPropertySupport</code> than determines the
 * above defined integer vectors <code>M(z|B)</code> and <code>N(z|B)</code>. It further
 * computes the transformed argument
 * <p align=center>
 * <code>
 * T(z|B) = z + 2&pi;iN(z|B) + BM(z|B)
 * </code>
 * <p>
 * and the exponential transformation factor
 * <p align=center>
 * <code>
 * f(z|B) = &frac12; M(z|B)<sup>t</sup>BM(z|B) + M(z|B)<sup>t</sup>z
 * </code>
 * <p>
 * which can be accessed by calling {@link #getTransformedZ()} and {@link #getFactor()}.
 * <p>
 * This class can also be used to check whether a vector <code>z&isin;C<sup>g</sup></code>
 * is a grid point of
 * <p align=center>
 * <code>
 * {2&pi;iN+BM|M,N&isin;Z<sup>g</sup>}.
 * </code>
 * <p>
 * In that case <code>T(z|B)</code> will obviously be zero.
 * @see Theta
 * @author Markus Schmies
 */
final public class TransformPropertySupport implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    final ComplexMatrix B = new ComplexMatrix();

    RealMatrix reBInv;

    final RealMatrix reB = new RealMatrix ();
    final RealMatrix imB = new RealMatrix ();

    final ComplexVector transfromedZ = new ComplexVector();

    final RealVector reBInvReZ      = new RealVector ();
    final RealVector M              = new RealVector ();
    final RealVector N              = new RealVector ();
    final RealVector tmp            = new RealVector ();

    final RealVector zRe            = new RealVector ();
    final RealVector zIm            = new RealVector ();

    final RealVector transformedZRe = new RealVector ();
    final RealVector transformedZIm = new RealVector ();

    final Complex factor = new Complex();

    final double minus2Pi = -2 * Math.PI;
    final double plus2Pi  =  2 * Math.PI;

    int dim;

    TransformPropertySupport(final int dim) {
        setDim(dim);
    }

    /**
     * Creates instance with presrcibed period matrix.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public TransformPropertySupport( final ComplexMatrix periodMatrix ) {
        setPeriodMatrix( periodMatrix );
    }

    void setDim( final int dim ) {
        
	if( this.dim != dim ) {

            transfromedZ.newSize(dim);

            M.newSize(dim);
            N.newSize(dim);

            tmp.newSize(dim);

            zRe.newSize(dim);
            zIm.newSize(dim);
            reB.newSize(dim);
            imB.newSize(dim);

            transformedZRe.newSize(dim);
            transformedZIm.newSize(dim);

            this.dim = dim;
        }
    }

    /**
     * Returns period matrix.
     */
    public ComplexMatrix getPeriodMatrix() {
        return new ComplexMatrix(B);
    }

    /**
     * Sets the period matrix.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public final void setPeriodMatrix( final ComplexMatrix periodMatrix ) {
        setDim(periodMatrix.getNumCols());

        B.assign(periodMatrix);

        B.getRe(reB);
        B.getIm(imB);

        reBInv = reB.invert();
    }

    /**
     * Returns the argument vector which was transformed.
     */
    public final ComplexVector getZ() {
        return new ComplexVector(zRe.re, zIm.re);
    }

    /**
     * Sets the argument vector which will be transformed.
     * @param z argument vector
     */
    public final void setZ(final ComplexVector z) {
        zRe.re = z.re;
        zIm.re = z.im;

        transformedZRe.re = transfromedZ.re;
        transformedZIm.re = transfromedZ.im;

        reBInvReZ.assignTimes(reBInv, zRe);

        M.assignRound(reBInvReZ);
        M.assignNeg();

        N.assignTimes(imB, M);
        N.assignPlus(zIm);
        N.assignTimes(1 / minus2Pi);
        N.assignRound();

        transformedZRe.assignTimes(reB, M);
        transformedZRe.assignPlus(zRe);
        transformedZIm.assignTimes(imB, M);
        transformedZIm.assignPlus(zIm);

        tmp.assignTimes(N, plus2Pi);
        transformedZIm.assignPlus(tmp);

        tmp.assignTimes(reB, M);
        factor.re = tmp.dot(M) / 2 + zRe.dot(M);

        tmp.assignTimes(imB, M);
        factor.im = tmp.dot(M) / 2 + zIm.dot(M);
    }

    /**
     * Returns the transformed argument vector.
     * @return <code> T(z|B) = z + 2&pi;iN(z|B) + BM(z|B) </code>
     */
    public final ComplexVector getTransformedZ() {
        return new ComplexVector(transfromedZ);
    }

    /**
     * Returns the exponential transformation factor.
     * @return <code> f(z|B) = &frac12; M(z|B)<sup>t</sup>BM(z|B) + M(z|B)<sup>t</sup>z </code>
     */
    public final Complex getFactor() {
        return new Complex(factor);
    }

    /**
     * Returns the integer vector <code>M(z|B)</code>.
     * @return M(z|B) = - [ re(B)<sup>-1</sup> re(z) ]
     */
    public final IntegerVector getM() {
        return IntegerVector.round(M);
    }

    /**
     * Returns the integer vector <code>N(z|B)</code>.
     * @return N(z|B) = - [ &frac12 ( im(z) + im(B)M(z|B )/&pi; ]
     */
    public final IntegerVector getN() {
        return IntegerVector.round(N);
    }
}
