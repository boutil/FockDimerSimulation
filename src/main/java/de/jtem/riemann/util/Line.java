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

package de.jtem.riemann.util;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexValue;

public class Line implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    public Complex P, Q;

    protected static double EPS = 1e-12;

    private static Complex U = new Complex();
    private static Complex V = new Complex();
    private static Complex W = new Complex();

    public Line() {
    	P = new Complex();
    	Q = new Complex();
    }

    public Line( ComplexValue aP, ComplexValue aQ ) {

	P = new Complex(aP);
	Q = new Complex(aQ);
    }

    public Line set( Complex aP, Complex aQ ) {
	P.assign(aP);
	Q.assign(aQ);
	return this;
    }
    public Line set( ComplexValue aP, ComplexValue aQ ) {
    	P.assign(aP);
    	Q.assign(aQ);
    	return this;
    }
    
    public boolean intersects( Line L ) {
	if( intersect( this, L, null, null ) != 0 )
	    return true;
	return false;
    }

    public boolean intersects( Line L, Complex I ) {
	if( intersect( this, L, I, null ) != 0 )
	    return true;
	return false;
    }

    public final double intersects( Line L2, Complex I, Complex coords ) {

	return intersect( this, L2, I, coords );
    }

    public static final double intersect( Line L1, Line L2, Complex I, Complex coords ) {
	return intersect( L1.P, L1.Q, L2.P, L2.Q, I, coords );
    }

    public static final double intersect( Complex l1P, Complex l1Q,
					  Complex l2P, Complex l2Q, Complex I, Complex coords ) {

	U.assignMinus( l2P, l1P );
	V.assignMinus( l1Q, l1P );
	W.assignMinus( l2P, l2Q );

	double det = V.re*W.im - V.im*W.re;

	if( det < EPS && det > - EPS )
	    return 0;

	if( coords != null  || I != null ) {

	    double coords0 = ( W.im*U.re - W.re*U.im ) / det;
	    double coords1 = (-V.im*U.re + V.re*U.im ) / det;

	    if( I != null ) {

		I.re = coords0 * l1Q.re + (1 - coords0) * l1P.re;
		I.im = coords0 * l1Q.im + (1 - coords0) * l1P.im;
	    }

	    if( coords != null ) {
		coords.re = coords0;
		coords.im = coords1;
	    }
	}

	return -det;
    }
}
