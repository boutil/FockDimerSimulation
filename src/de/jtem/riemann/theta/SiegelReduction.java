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
import de.jtem.blas.IntegerMatrix;
import de.jtem.blas.RealMatrix;
import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.algebra.linear.decompose.Cholesky;
import de.jtem.numericalMethods.geometry.latticeReduction.LLL;

/**
 * Provides Siegel`s reduction algorithm.
 *
 * The modular transformation property of Riemann theta function
 * enables us to choose a period matrix from a hole family given by all
 * <p align=center>
 * <code>
 * <table>
 * <tr>
 * <td> &sigma;(B) = 2&pi;i(aB+2&pi;ib)(cB+2&pi;id)<sup>-1</sup>,
 * </td>
 * <td> where &sigma; =
 * </td>
 * <td>
 * <table>
 * <tr> <td>/ a b \</td><td> </td></tr>
 * <tr> <td>\ c d /</td><td>,</td></tr>
 * </table>
 * </td>
 * </tr>
 * </table>
 * </code>
 * <p>
 * is any element of the {@link ModularGroup}. This choice is one-to-one
 * related to that of the homology basis for a Riemann surface.
 * <p>
 * The question is now: which one is best? E.g. for which element of the
 * modular group (resp. homology basis) converges the Fourier series of the
 * Riemann theta function most rapidly? Siegel's Reduction algorithm is
 * an answer to this, although it is not optimal.
 * One major problem of its implementation is that it requires to perform
 * a lattice reduction, which is know to be NP hard.
 * However, an approximate algorithm due to Lenstra, Lenstra and Lovasz,
 * the so called {@link de.jtem.numericalMethods.geometry.latticeReduction.LLL} algorithm
 * gives satisfactory and fast results.
 * <p>
 * @see ModularPropertySupport
 * @see ModularTransformation
 * @see de.jtem.numericalMethods.geometry.latticeReduction.LLL
 * @see Theta
 * @author Markus Schmies
 */
public class SiegelReduction implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    final static Complex PI2i = new Complex(0, 2 * Math.PI);

    final static double EPS = 1e-13;

    int dim;

    ModularTransformation modular, specialModularTransformation, tmpModular;

    ComplexMatrix PM, tPM;

    RealMatrix tX;
    RealMatrix tY;
    RealMatrix tL;
    RealMatrix tLInv;

    IntegerMatrix A;
    IntegerMatrix B;
    IntegerMatrix AInvTr;
    IntegerMatrix tmpIntegerMatrix;

    RealMatrix doubleA      = new RealMatrix ();
    RealMatrix doubleAInvTr = new RealMatrix ();

    ComplexMatrix tmpComplexMatrix = new ComplexMatrix();
    ComplexMatrix cD, cN, cDInv; // temp matrices are use in apply

    final Complex tPM00      = new Complex();
    final Complex tmpComplex = new Complex();

    IntegerMatrix U;

    boolean forNegativeDefinitePeriodMatrices = true;

    SiegelReduction(int dim) {
        modular = new ModularTransformation(dim);

        PM = new ComplexMatrix(dim);
        tPM = new ComplexMatrix(dim);

        PM.setDiagonal(2 * Math.PI, 0);
        tPM.setDiagonal(2 * Math.PI, 0);

        setDim(dim);
    }

    /**
     * Creates instance with presribed period matrix and
     * performs Siegel's reduction algorithm.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     */
    public SiegelReduction( final ComplexMatrix periodMatrix ) {
        this( periodMatrix.getNumRows() );
        setPeriodMatrix( periodMatrix );
    }

    private void setDim( final int dim) {

        if (this.dim != dim) {
            this.dim = dim;

            tX = new RealMatrix (dim);
            tY = new RealMatrix (dim);
            tL = new RealMatrix (dim);

            tLInv = new RealMatrix (dim);

            A            = new IntegerMatrix (dim);
            AInvTr       = new IntegerMatrix (dim);
            B            = new IntegerMatrix (dim);
            tmpIntegerMatrix = new IntegerMatrix (dim);
            U            = new IntegerMatrix (dim);

            tmpModular = new ModularTransformation(dim);

            specialModularTransformation = getSpecialModularTransformation(dim);

            cD    = new ComplexMatrix();
            cN    = new ComplexMatrix();
            cDInv = new ComplexMatrix();
        }
    }

    /**
     * Returns <code>&sigma;</code>, the  element of the modular group,
     * determing <code>&sigma;(B)</code>, the result of Siegel's reduction algorithm.
     * <p align=center>
     * <code>
     * <table>
     * <tr>
     * <td> &sigma;(B) = 2&pi;i(aB+2&pi;ib)(cB+2&pi;id)<sup>-1</sup>,
     * </td>
     * <td> where &sigma; =
     * </td>
     * <td>
     * <table>
     * <tr> <td>/ A B \</td><td> </td></tr>
     * <tr> <td>\ C D /</td><td>,</td></tr>
     * </table>
     * </td>
     * </tr>
     * </table>
     * </code>
     * @return <code>&sigma;</code>, the  element of the modular group realizing the result
     * of the reduction
     */
    public ModularTransformation getModularTransformation() {
        return new ModularTransformation(modular);
    }

    /**
     * Returns period matrix which was reduced.
     */
    public ComplexMatrix getPeriodMatrix() {
        return new ComplexMatrix(PM);
    }

    /**
     * Sets period matrix which will be reduces.
     */
    public void setPeriodMatrix( final ComplexMatrix aPeriodMatrix ) {
        PM.assign(aPeriodMatrix);
        compute();
    }

    /**
     * Returns the reduced period matrix <code>&sigma;(B)</code>.
     * <p align=center>
     * <code>
     * <table>
     * <tr>
     * <td> &sigma;(B) = 2&pi;i(aB+2&pi;ib)(cB+2&pi;id)<sup>-1</sup>,
     * </td>
     * <td> where &sigma; =
     * </td>
     * <td>
     * <table>
     * <tr> <td>/ A B \</td><td> </td></tr>
     * <tr> <td>\ C D /</td><td>,</td></tr>
     * </table>
     * </td>
     * </tr>
     * </table>
     * </code>
     * @return <code>&sigma;(B)</code>, the reduced period matrix.
     */
    public ComplexMatrix getReducedPeriodMatrix() {
        return new ComplexMatrix(tPM);
    }

    boolean shiftImaginaryPart() {
        boolean shifted = false;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                int shift = -(int)Math.floor(tPM.im[i][j] / 2 / Math.PI + 0.5);
                if (shift != 0) {
                    shifted = true;
                    B.re[i][j] = shift;
                } else {
                    B.re[i][j] = 0;
                }
            }
        if (!shifted)
            return false;
        tmpComplexMatrix.assignTimes(B, PI2i);
        tPM.assignPlus(tmpComplexMatrix);
        modular.applyGenerator(B);

        return true;
    }

    boolean reduceLattice() {
        tPM.getRe(tX);
        Cholesky.decompose(tX.re, tL.re);
        tLInv.assignInvert(tL);
        LLL.reduce(tL.re, U.re); // the lll impl. takes an array of vectors as a basis, not a matrix
        if (U.isId())
            return false;
        doubleA.assignTimes(tL, tLInv);
        A.assignRound(doubleA);
        doubleAInvTr.assignInvert(doubleA);
        doubleAInvTr.assignTranspose();
        AInvTr.assignRound(doubleAInvTr);
        modular.applyGenerator(A, AInvTr);
        tmpComplexMatrix.assignTimes(A, tPM);
        A.assignTranspose();
        tPM.assignTimes(tmpComplexMatrix, A);
        return true;
    }

    boolean applySpecialTransformation() {
        tPM.get(0, 0, tPM00);

        if (tPM00.absSqr() + EPS > 4 * Math.PI * Math.PI)
            return false;

        apply(specialModularTransformation);
        return true;
    }

    void apply(ModularTransformation M) {
        tmpModular.assign(modular); modular.assignTimes(M, tmpModular);
        if (!modular.respectsSymplecticStructure()) {
            System.out.println(modular);
            throw new RuntimeException(" does not ");
        }
        computeReducedPeriodMatrix();
    }

    void computeReducedPeriodMatrix() {
        tPM.assign(PM);
        applyToReducedPeriodMatrix(modular);
    }

    void applyToReducedPeriodMatrix( final ModularTransformation M ) {
        cN.assignTimes(M.a, tPM);
        cN.assignDivide(PI2i);
        cN.assignPlus(M.b);
        cD.assignTimes(M.c, tPM);
        cD.assignDivide(PI2i);
        cD.assignPlus(M.d);
        cDInv.assignInvert(cD);
        tPM.assignTimes(cN, cDInv);
        tPM.assignTimes(PI2i);
    }

    void checkForNegativeDefinitePeriodMatrices() {
        if (forNegativeDefinitePeriodMatrices) {
            PM.assignTimes(-1);
            modular.b.assignTimes(-1);
            modular.c.assignTimes(-1);
        }
    }

    void compute() {
        checkForNegativeDefinitePeriodMatrices();
        computeReducedPeriodMatrix();
        try {
            boolean wasChanged = true;
            while (wasChanged) {
                wasChanged = reduceLattice();
                wasChanged = applySpecialTransformation() || wasChanged;
                wasChanged = shiftImaginaryPart() || wasChanged;
            }
        } finally {
            checkForNegativeDefinitePeriodMatrices();
        }
        computeReducedPeriodMatrix();
    }

    static ModularTransformation getSpecialModularTransformation( final int dim ) {
        IntegerMatrix a = new IntegerMatrix (dim);
        IntegerMatrix b = new IntegerMatrix (dim);
        IntegerMatrix c = new IntegerMatrix (dim);
        IntegerMatrix d = a;

        a.assignId();
        a.set(0, 0, 0);
        b.set(0, 0, -1);
        c.set(0, 0, 1);

        return new ModularTransformation(a, b, c, d);
    }
}
