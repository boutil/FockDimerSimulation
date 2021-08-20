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

package de.jtem.riemann.schottky;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;

public class SchottkyDomainSampler implements Serializable{

    private static final long serialVersionUID = 1L;

    public static double [] getBound( Schottky schottky )  {

	double minX, maxX, minY, maxY;

	minX = maxX = schottky.center[0][0].re;
	minY = maxY = schottky.center[0][0].im;

	for( int i=0; i<schottky.numGenerators; i++ )
	    for( int j=0; j<2;j++ ) {

		if( schottky.center[i][j].re - schottky.radius[i] < minX )
		    minX = schottky.center[i][j].re - schottky.radius[i];

		if( schottky.center[i][j].im - schottky.radius[i] < minY )
		    minY = schottky.center[i][j].im - schottky.radius[i];

		if( schottky.center[i][j].re + schottky.radius[i] > maxX )
		    maxX = schottky.center[i][j].re + schottky.radius[i];

		if( schottky.center[i][j].im + schottky.radius[i] > maxY )
		    maxY = schottky.center[i][j].im + schottky.radius[i];
	    }

	return new double[] { minX, minY, maxX, maxY };
    }


    public static double [] getRaster( Schottky schottky,
				       int numRows, int numCols, double xDist, double yDist ) {

	double [] bound = getBound( schottky );

	bound[0] -= xDist;
	bound[1] -= yDist;
	bound[2] += xDist;
	bound[3] += yDist;

	return getRaster( schottky, numRows, numCols, bound );
    }

    public static double [] getRaster( Schottky schottky,
				       int numRows, int numCols, double [] bound ) {

	Complex z = new Complex();

	double [] xy = new double[2*(numRows+1)*(numCols+1)];

	double deltaX = (bound[2]-bound[0]) / numCols;
	double deltaY = (bound[3]-bound[1]) / numRows;

	double row = bound[1];

	int k=0;

	for( int i=0; i<=numRows; i++, row += deltaY ) {

	    double col = bound[0];

	    for( int j=0; j<=numCols; j++, col += deltaX ) {

		z.assign( col, row );

		if( schottky.dist( z ) == 0 ) {
		    xy[k++] = col;
		    xy[k++] = row;
		}
	    }
	}

	double [] rasterXY = new double[k];

	System.arraycopy( xy, 0, rasterXY, 0, k );

	return rasterXY;
    }

    public static double [] getCover( Schottky schottky, double dense, double xDist, double yDist ) {

	double [] bound = getBound( schottky );

	bound[0] -= xDist;
	bound[1] -= yDist;
	bound[2] += xDist;
	bound[3] += yDist;

	return getCover( schottky, dense, bound );
    }

    public static double [] getCover( Schottky schottky, double dense, double [] bound ) {
	return getCover( schottky, dense, dense, bound );
    }

    public static double [] getCover( Schottky schottky,
				      double dense, double ringDense, double [] bound ) {

	int numRows = (int)(1 / dense);
	int numCols = (int)(1 / dense);

	double [] raster = getRaster( schottky, numRows, numCols, bound );

	int ringLength = 0;

	for( int i=0; i<schottky.numGenerators; i++ )
	    ringLength += 2 * Math.max( 4, 2*(int)(1 / ringDense) );

	double [] cover = new double [ raster.length + 2 * ringLength ];

	System.arraycopy( raster, 0, cover, 0, raster.length );

	int k = raster.length;

	for( int i=0; i<schottky.numGenerators; i++ ) {

	    double numPointsOnRing = Math.max( 4, 2*(int)(1 / ringDense) );

	    for( int j=0; j<2; j++ ) {

		for( int l=0; l<numPointsOnRing; l++ ) {
		    cover[k  ] = schottky.center[i][j].re
			+ schottky.radius[i] * Math.cos( 2 * l/numPointsOnRing * Math.PI ) * 1.0001;
		    cover[k+1] = schottky.center[i][j].im
			+ schottky.radius[i] * Math.sin( 2 * l/numPointsOnRing * Math.PI ) * 1.0001;

		    if( cover[k  ] >= bound[0] && cover[k  ]<=bound[2] &&
			cover[k+1] >= bound[1] && cover[k+1]<=bound[3]   )
			k+=2;
		}
	    }
	}

	double [] finalCover = new double[k];

	System.arraycopy( cover, 0, finalCover, 0, k );

	return finalCover;
    }

}





