// JTEM - Java Tools for Experimental Mathematics
// Copyright (C) 2001 JEM-Group
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

package de.jtem.riemann.constrainedComplex;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.SurfacePoint;

public class Real extends SingleConstrain {
    
    private static final long serialVersionUID = 1L;

    public Real( SurfacePoint aA ) {

	A = aA;

	A.setConstrain( this );
    }

    Complex real = new Complex();

    public void projectedMove( ConstrainedComplex aConstrainedComplex, Complex newCoords ) {
	real.re = newCoords.re;

	move( aConstrainedComplex, real );
    }

    public void move( ConstrainedComplex aConstrainedComplex, Complex newCoords ) {

	if( Math.abs( newCoords.im ) > EPS || aConstrainedComplex != A )
	    throw new IllegalCoordinateException( aConstrainedComplex );

	real.re = newCoords.re;

	A.move( real );
    }

    public int getNumOfParameters( ConstrainedComplex aConstrainedComplex ) {
	return aConstrainedComplex == A ? 1 : 0;
    }
    
    public void setByParameter( double [] p , int offset, ConstrainedComplex aConstrainedComplex ) {
	if( aConstrainedComplex == A ) {
	    real.re = p[offset];
	    move( A, real );
	}
    }

    public void getValue( double[] p, int offset, ConstrainedComplex aConstrainedComplex ) {
	if( aConstrainedComplex == A ) {
	    p[offset  ] = A.coords.re;
	}
    }

}












