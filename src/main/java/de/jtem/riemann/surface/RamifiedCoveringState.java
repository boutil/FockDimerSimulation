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

package de.jtem.riemann.surface;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.homologie.TransformState;

public class RamifiedCoveringState implements java.io.Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    private RamifiedCovering surface;

    private Complex [] initSurfacePointCoords;

    private TransformState transformState;

    RamifiedCoveringState( RamifiedCovering aSurface ) {
	surface = aSurface;

	initSurfacePointCoords = new Complex     [ surface.numOfSurfacePoints ];

	for( int i=0; i < surface.numOfSurfacePoints; i++ )
	    initSurfacePointCoords[i] = new Complex();

	save( surface );
    }

    void save( RamifiedCovering sender ) {

	if( sender != surface )
	    throw new IllegalArgumentException( "surface state does not belong to sender" );

	if( surface.getTransform() != null ) {
	    if( transformState == null )
		transformState = surface.getTransform().getState();
	    else
		surface.getTransform().getState( transformState );
	} else {
	    transformState = null;
	}

	for( int i=0; i < surface.numOfSurfacePoints; i++ ) {
	    initSurfacePointCoords[i].assign( surface.initSurfacePoint[i] );
	}
    }

    void load( RamifiedCovering sender ) {

	if( sender != surface )
	    throw new IllegalArgumentException( "surface state does not belong to sender" );

	if( surface.getTransform() != null )
	    surface.getTransform().setState( transformState );

	for( int i=0; i < surface.numOfSurfacePoints; i++ ) {
	    surface.initSurfacePoint[i].setNewCoords( initSurfacePointCoords[i] );
	}

	java.util.Arrays.sort( surface.surfacePoint );

	surface.listOfPasses.clear();

    }
}













