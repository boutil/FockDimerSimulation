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

public final class LineSeg extends Line implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;
 
    public LineSeg() {
    }

    public LineSeg( ComplexValue aP, ComplexValue aQ ) {
    	super( aP, aQ );
    }

    public boolean intersects( LineSeg L ) {
	if( intersect( this, L, null, null ) != 0 )
	    return true;
	return false;
    }


    public boolean intersects( LineSeg L, Complex I ) {
	if( intersect( this, L, I, null ) != 0 )
	    return true;
	return false;
    }

    public boolean intersects( Line L, Complex I ) {
	if( intersect( L, this, I, null ) != 0 )
	    return true;
	return false;
    }

    public boolean intersects( Line L ) {
	if( intersect( L, this, null, null ) != 0 )
	    return true;
	return false;
    }

    public final double intersects( LineSeg L2, Complex I, Complex coords ) {

	return intersect( this, L2, I, coords );
    }

    private static Complex tmpCoords = new Complex();
 
    static double intersect( LineSeg L1, LineSeg L2, Complex I, Complex coords ) {
	
	if( coords == null )
	    coords = new Complex();

	if( I == null )
	    I = new Complex();

	double det = Line.intersect( L1, L2, I, coords );
	
	if( det != 0
	    && coords.re < 1+EPS && coords.re > -EPS
	    && coords.im < 1+EPS && coords.im > -EPS ) return det;

	return 0;
    }

    static double intersect( Line L1, LineSeg L2, Complex I, Complex coords ) {
	
	if( coords == null )
	    coords = new Complex();

	double det = Line.intersect( L1, L2, I, coords );
	
	if( det != 0
	    && coords.im < 1+EPS && coords.im > -EPS ) return det;

	return 0;
    }

}
