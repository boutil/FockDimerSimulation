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

import java.io.Serializable;

import de.jtem.mfc.field.Complex;

public abstract class Constrain implements Serializable, Cloneable {
    
    private final static long serialVersionUID = -5243844056110607340L;

    final double EPS = 1e-12;

    /**
     * 
     * @param aConstrainedComplex
     * @param newCoords
     */
    public void projectedMove( ConstrainedComplex aConstrainedComplex, Complex newCoords ) {
    	move( aConstrainedComplex, newCoords );
    }

    public abstract void move( ConstrainedComplex aConstrainedComplex, Complex newCoords )
		throws IllegalCoordinateException; 
    
    public abstract boolean effects( ConstrainedComplex aConstrainedComplex );

    public abstract void deconstrain();


    public abstract int getNumOfParameters( ConstrainedComplex aPoint );

    public abstract void setByParameter( double [] p, int offset, ConstrainedComplex aConstrainedComplex );
    public abstract void getValue(       double [] p, int offset, ConstrainedComplex aConstrainedComplex );
}












