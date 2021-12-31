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

import java.awt.Color;
import java.util.Vector;

import de.jtem.moebiusViewer.Attributes;
import de.jtem.moebiusViewer.MoebiusGraphics;
import de.jtem.moebiusViewer.shape.AbstractShape;
import de.jtem.moebiusViewer.shape.Point;
import de.jtem.moebiusViewer.shape.Rectangle;


public class RamifiedCoveringShape extends AbstractShape {

    private static final long serialVersionUID = 1L;

    public Point [] point;

    protected RamifiedCovering surface;

    /**
       * Get the value of surface.
       * @return Value of surface.
       */
    public RamifiedCovering getSurface() {return surface;}

    /**
       * Set the value of surface.
       * @param v  Value to assign to surface.
       */
    public void setSurface(RamifiedCovering  v) {
	set( v );
    }

    public void set( RamifiedCovering v ) {

	if( surface != v ) {

	    surface = v;

	    point = new Point[surface.numOfSurfacePoints + 1];

	    point[0] = new Point();

	    point[0].setColor( Color.white );

	    int j = 1;

	    for( int i=0; i<surface.numOfBranchPoints; i++, j++ ) {

		point[j] = new Point();
		point[j].setColor( Color.black );
	    }

	    for( int i=0; i<surface.numOfSingularPoints; i++, j++ ) {

		point[j] = new Point();
		point[j].setColor( Color.red );
	    }

	    for( int i=0; i<surface.numOfDistinguishedPoints; i++, j++ ) {

		point[j] = new Point();
		point[j].setColor( Color.blue );
	    }
	}

	compute();

	firePropertyChange( "surface" );
    }


    public void compute() {

	point[0].setCoords( surface.origin.re, surface.origin.im );

	int j = 1;

	for( int i=0; i<surface.numOfBranchPoints; i++, j++ ) {
	    point[j].setCoords( surface.branchPoint[i] );
	}

	for( int i=0; i<surface.numOfSingularPoints; i++, j++ ) {

	    point[j].setCoords( surface.singularPoint[i] );
	}

	for( int i=0; i<surface.numOfDistinguishedPoints; i++, j++ ) {

	    point[j].setCoords( surface.distinguishedPoint[i] );
	}
    }

    public RamifiedCoveringShape() {
    }

    public RamifiedCoveringShape( RamifiedCovering aSurface ) {
	setSurface( aSurface );
    }

    public void draw (MoebiusGraphics context) {

	compute();

	context.circle( 0, 0, 1 );

	if( surface == null )
	    return;

	for( int i=0; i<point.length; i++ )
	    point[i].draw( context );

	context.setVerticalTextLayout(   Attributes.TOP );
	context.setHorizontalTextLayout( Attributes.RIGHT );

	context.setColor( Color.black );

	for( int i=0; i<surface.numOfBranchPoints; i++ )

	  context.text( surface.branchPoint[i].getRe(),
	  	  surface.branchPoint[i].getIm(), " " + i );

	context.setColor( Color.red );

	for( int i=0; i<surface.numOfSingularPoints; i++ )
	    context.text( surface.singularPoint[i].getRe(),
			  surface.singularPoint[i].getIm(), " " + i);

	context.setVerticalTextLayout( Attributes.BOTTOM );

	context.setColor( Color.blue );

	for( int i=0; i<surface.numOfDistinguishedPoints; i++ )

	  context.text( surface.distinguishedPoint[i].getRe(),
	  	  surface.distinguishedPoint[i].getIm(), i + " " );

    }

    public Rectangle getBound() {

	if( point.length < 1 )
	    return null;

	double maxX = 1;
	double maxY = 1;

	for( int i = 0; i<point.length; i++ ) {
	    double aX = Math.abs( point[i].getX() );
	    double aY = Math.abs( point[i].getY() );

	    maxX = maxX<aX ? aX : maxX;
	    maxY = maxY<aY ? aY : maxY;
	}

	return new Rectangle( -maxX, -maxY, 2*maxX, 2*maxY );
    }

    public Vector getTools() {

	Vector toolList = new Vector();

	toolList.addElement( new DragSurfacePointsTool() );

	return toolList;
    }

}












