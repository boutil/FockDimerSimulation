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

public class RealSymmetry extends ToupleConstrain {

    private static final long serialVersionUID = 1L;

    public RealSymmetry( SurfacePoint aA, SurfacePoint aB ) {
	super( aA, aB );
    }

    public void move( ConstrainedComplex aConstrainedComplex, Complex newCoords ) {

//	double absSqrOfNewCoords = newCoords.absSqr();

	if( newCoords.im == 0 )
	    throw new IllegalCoordinateException( aConstrainedComplex );

	if( newCoords.im > 0 && A.coords.im > 0 || newCoords.im < 0 && A.coords.im < 0 ) {
	    newCoordsOfA.assign( newCoords );
	    newCoordsOfB.assignConjugate( newCoords );
	} else {
	    newCoordsOfB.assign( newCoords );
	    newCoordsOfA.assignConjugate( newCoords );
	}

	A.move( newCoordsOfA );
	B.move( newCoordsOfB );
    }
}












