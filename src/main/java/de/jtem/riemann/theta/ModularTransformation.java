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

import de.jtem.blas.IntegerMatrix;

/**
 * An instance of this class represents an element of the modular group.
 *
 * The modular group consists of square integer matrices
 * <p align=center>
 * <code>
 * <table>
 * <tr> 
 * <td> &sigma; =
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
 * <p>
 * with <code>A,B,C,D&isin;Z<sup>g<small>&times;</small>g</sup></code>, 
 * that leave the symplectic structure invariant, e.g.:
 * <p align=center>
 * <code>
 * <table>
 * <tr> 
 * <td>
 * &sigma;<sup>t</sup>
 * </td>
 * <td>
 * <table>
 * <tr> <td><pre>/ 0 1 \</pre></td></tr>
 * <tr> <td><pre>\-1 0 /</pre></td></tr>
 * </table>
 * </td>
 * <td>
 * &sigma; = 
 * </td>
 * <td>
 * <table>
 * <tr> <td><pre>/ 0 1 \</pre)</td></tr>
 * <tr> <td><pre>\-1 0 /</pre></td></tr>
 * </table>
 * </td>
 * </tr>
 * </table>
 * </code>
 *
 * @see ModularPropertySupport
 * @see SiegelReduction
 * @author Markus Schmies
 */
public class ModularTransformation implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    int dim;

    IntegerMatrix a;
    IntegerMatrix b;
    IntegerMatrix c;
    IntegerMatrix d;

    private IntegerMatrix at;
    private IntegerMatrix bt;
    private IntegerMatrix ct;
    private IntegerMatrix dt;
 // transposed of at, bt, ct, dt used in getDefect

    private IntegerMatrix p;
    private IntegerMatrix q;
    private IntegerMatrix e;
        // p and q are temp matrices, e is the identity

    /**
     * Creates identity in molular group of degree <code>2g</code>.
     */
    public ModularTransformation( final int g ) {
        a = IntegerMatrix.id(g);
        d = IntegerMatrix.id(g);

        b = new IntegerMatrix (g);
        c = new IntegerMatrix (g);

        setDim(g);
    }

    /**
     * Creates a copy of <code>T</code>.
     */
    public ModularTransformation( final ModularTransformation T ) {
        a = new IntegerMatrix (T.a);
        b = new IntegerMatrix (T.b);
        c = new IntegerMatrix (T.c);
        d = new IntegerMatrix (T.d);

        setDim(T.dim);
    }

    /**
     * Creates element with prescribed sub-matrices <code>A,B,C</code> and, <code>D</code>.
     * @param A square integer matrix of degree <code>dim</code>
     * @param B square integer matrix of degree <code>dim</code>
     * @param C square integer matrix of degree <code>dim</code>
     * @param D square integer matrix of degree <code>dim</code>
     */
    public ModularTransformation( final IntegerMatrix A, final IntegerMatrix B, final IntegerMatrix C, final IntegerMatrix D ) {
        A.checkSquared();
        B.checkShape(A);
        C.checkShape(A);
        D.checkShape(A);

        a = new IntegerMatrix (A);
        b = new IntegerMatrix (B);
        c = new IntegerMatrix (C);
        d = new IntegerMatrix (D);

        setDim(A.getNumRows());

        if (!respectsSymplecticStructure())
            throw new IllegalArgumentException("data do not respect symplectic structure");
    }

    /**
     * Queries sub-matrix a.
     * @return sub-matrix a
     */
    public IntegerMatrix getA() {
    	return new IntegerMatrix(a);
    }
    /**
     * Queries sub-matrix b.
     * @return sub-matrix b
     */
    public IntegerMatrix getB() {
    	return new IntegerMatrix(b);
    }
    /**
     * Queries sub-matrix c.
     * @return sub-matrix c
     */
    public IntegerMatrix getC() {
    	return new IntegerMatrix(c);
    }
    /**
     * Queries sub-matrix d.
     * @return sub-matrix d
     */
    public IntegerMatrix getD() {
    	return new IntegerMatrix(d);
    }
    
    /**
     * Returns degree of modular group.
     */
    public int getDegree() {
        return 2 * dim;
    }

    private void setDim(int aDim) {
        if (aDim != dim) {

            dim = aDim;

            at = new IntegerMatrix (aDim);
            bt = new IntegerMatrix (aDim);
            ct = new IntegerMatrix (aDim);
            dt = new IntegerMatrix (aDim);

            p = new IntegerMatrix (dim);
            q = new IntegerMatrix (dim);

            e = IntegerMatrix.id(dim);
        }
    }

    /**
     * Assigns <code>this</code> with identity in modular group.
     */
    public void assignId() {
        a.assignId();
        b.assignZero();
        c.assignZero();
        d.assignId();
    }
    
    /**
     * Tests if <code>this</code> is the identity in modular group.
     */
    public boolean isId() {
        return a.isId() && b.isZero() && c.isZero() && d.isId();
    }

    /**
     * Assigns <code>this</code> with <code>T</code>.
     */
    public void assign( final ModularTransformation T ) {
        a.assign(T.a);
        b.assign(T.b);
        c.assign(T.c);
        d.assign(T.d);
        setDim(T.dim);
    }

    /**
     * Assigns <code>this</code> with inverse of <code>T</code>.
     */
    public void assignInvert( final ModularTransformation T ) {

	T.computeTransposed();

        a.assign   (T.dt);
        b.assignNeg(T.bt);
        c.assignNeg(T.ct);
        d.assign   (T.at);
	
        setDim(T.dim);
    }

    /**
     * Returns inverse of <code>this</code>.
     */
    public ModularTransformation invert() {
	ModularTransformation T = new ModularTransformation( dim );
	T.assignInvert( this );
	return T;
    }
    
    void computeTransposed() {

	at.assignTranspose(a);
	bt.assignTranspose(b);
        ct.assignTranspose(c);
	dt.assignTranspose(d);
    }

    boolean respectsSymplecticStructure() {
        return getDefect() < 1e-14;
    }

    double getDefect() {
        int defect = 0;

        bt.assignTranspose(b); p.assignTimes(a, bt);
        at.assignTranspose(a); q.assignTimes(b, at);

        p.assignMinus(q);

        defect += p.normSqr();

        dt.assignTranspose(d); p.assignTimes(a, dt);
        ct.assignTranspose(c); q.assignTimes(b, ct);

        p.assignMinus(q);
        p.assignMinus(e);

        defect += p.normSqr();

        p.assignTimes(c, bt);
        q.assignTimes(d, at);

        p.assignMinus(q);
        p.assignPlus(e);

        defect += p.normSqr();

        p.assignTimes(c, dt);
        q.assignTimes(d, ct);

        p.assignMinus(q);

        defect += p.normSqr();

        return defect;
    }

    /**
     * Returns product of <code>this</code> and <code>M</code>.
     * @param M element with the same degree as <code>this</code>
     */
    public ModularTransformation times( ModularTransformation M ) {
	ModularTransformation T = new ModularTransformation( dim );
	T.assignTimes( this, M );
	return T;
    }

    /**
     * Assigns <code>this</code> with product of itself and <code>M</code>.
     * @param M element with the same degree as <code>this</code>
     */
    public void assignTimes(ModularTransformation M) {
        assignTimes(this, M);
    }

    /**
     * Assigns <code>this</code> with product of <code>M</code> and <code>N</code>.
     * @param M element with the same degree as <code>N</code>
     * @param N element with the same degree as <code>M</code>
     */
    public void assignTimes(ModularTransformation M, ModularTransformation N) {
        if (M.dim != N.dim)
            throw new IllegalArgumentException("modular transformations have diffrent degree");

        setDim(M.dim);

        
    IntegerMatrix aM;
    IntegerMatrix bM;
    IntegerMatrix cM;
    IntegerMatrix dM;

        
    IntegerMatrix aN;
    IntegerMatrix bN;
    IntegerMatrix cN;
    IntegerMatrix dN;


        if (this == M) {
            aM = new IntegerMatrix (M.a); bM = new IntegerMatrix (M.b);
            cM = new IntegerMatrix (M.c); dM = new IntegerMatrix (M.d);
        } else {
            aM = M.a; bM = M.b;
            cM = M.c; dM = M.d;
        }
        if (M == N) {
            aN = aM; bN = bM;
            cN = cM; dN = dM;
        } else if (this == N) {
            aN = new IntegerMatrix (N.a); bN = new IntegerMatrix (N.b);
            cN = new IntegerMatrix (N.c); dN = new IntegerMatrix (N.d);
        } else {
            aN = N.a; bN = N.b;
            cN = N.c; dN = N.d;
        }

        p.assignTimes(aM, aN); q.assignTimes(bM, cN); a.assignPlus(p, q);
        p.assignTimes(aM, bN); q.assignTimes(bM, dN); b.assignPlus(p, q);
        p.assignTimes(cM, aN); q.assignTimes(dM, cN); c.assignPlus(p, q);
        p.assignTimes(cM, bN); q.assignTimes(dM, dN); d.assignPlus(p, q);
    }

    /**
     * Returns product of <code>this</code> and inverse of<code>M</code>.
     * @param M element with the same degree as <code>this</code>
     */
    public ModularTransformation divide( ModularTransformation M ) {
	ModularTransformation T = new ModularTransformation( dim );
	T.assignDivide( this, M );
	return T;
    }

    /**
     * Assigns <code>this</code> with product of itself and inverse of <code>M</code>.
     * @param M element with the same degree as <code>this</code>
     */
    public void assignDivide(ModularTransformation M) {
        assignDivide(this, M);
    }

    /**
     * Assigns <code>this</code> with product of <code>M</code> and inverse of <code>N</code>.
     * @param M element with the same degree as <code>N</code>
     * @param N element with the same degree as <code>M</code>
     */
    public void assignDivide(ModularTransformation M, ModularTransformation N) {
        if (M.dim != N.dim)
            throw new IllegalArgumentException("modular transformations have diffrent degree");
	
	this.assignInvert( N );
	this.assignTimes( M, this );
    }

    /**
     * Returns the for sub-matrices as string.
     */
    public String toString() {
        StringBuffer sb = new StringBuffer(300);

        sb.append("A=");
        sb.append(a);
        sb.append("\nB=");
        sb.append(b);
        sb.append("\nC=");
        sb.append(c);
        sb.append("\nD=");
        sb.append(d);

        return sb.toString();
    }

    /** applies generator of type 1 from the left which is b=c=0 and a transposed
	equals d inverse. Thus either D is null and d will be computed or
	D has already to equal the inverse of the transposed of A
	which is not checked. */
    void applyGenerator(IntegerMatrix A, IntegerMatrix D) {
        if (D == null) {
            D = A.invert();
            D.assignTranspose();
        }

        p.assign(a); a.assignTimes(A, p);
        p.assign(b); b.assignTimes(A, p);
        p.assign(c); c.assignTimes(D, p);
        p.assign(d); d.assignTimes(D, p);
    }

    /** applies generator of type 2 from the left which is a=d=id, c=0 and
	b symmetric. It is not check whether B is symmetric. */
    void applyGenerator(IntegerMatrix B) {
        q.assignTimes(B, c); a.assignPlus(q);
        q.assignTimes(B, d); b.assignPlus(q);
    }

    /** applies (from right) generator of type 3 which is the symplectic strucutre. */
    void applyGenerator() {
        IntegerMatrix t;

        t = c; c = a; a = t;
        t = d; d = b; b = t;

        a.assignTimes(-1);
        b.assignTimes(-1);
    }
}

