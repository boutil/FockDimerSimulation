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

import de.jtem.blas.RealMatrix;
import de.jtem.blas.RealVector;
import de.jtem.numericalMethods.algebra.linear.VectorOperations;

/**
 * @deprecated
 * @author schmies
 */
public final class LatticePointsForUniformApproximation_old 
    extends LatticePointsInEllipsoid_old implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    int volume;

    int [] size, product;

    RealVector origin;
    RealVector tmp, tmp_;

    boolean [] isSelected = new boolean[1000];

    RealVector[] delta;

    int deltaLength;

    double originOffset;

    boolean useSymmetrie = false;

    public LatticePointsForUniformApproximation_old( RealMatrix B ) {
	super( B );
    }

    public LatticePointsForUniformApproximation_old( RealMatrix B, double aRadius ) {	
	super( B, aRadius );
    }

    public LatticePointsForUniformApproximation_old( RealMatrix B, boolean useSymmetrie ) {	
	super( B );

	this.useSymmetrie = useSymmetrie;
    }

    public LatticePointsForUniformApproximation_old( RealMatrix B, double aRadius, boolean useSymmetrie ) {	
	super( B, aRadius );

	this.useSymmetrie = useSymmetrie;
    }



    private int numOfNewDelta;

    private void computeDelta() {

	deltaLength = (int)(Math.pow( 2, dim ) + 0.5 );

	delta = new RealVector[ deltaLength ];

	numOfNewDelta = 0;

	computeDelta( 0 );
    }

    
    private void computeDelta( int k ) {

	if( k != dim )  {

	    y.re[k] = 0; computeDelta( k+1 );
	    y.re[k] = 1; computeDelta( k+1 );

	} else delta[ numOfNewDelta++ ] = new RealVector ( y );	    
    }

    protected void setDim( int newDim ) { 

	if( newDim != dim ) {
	    
	    super.setDim( newDim );
	    
	    size    = new int[dim];

	    if( dim > 1 )
		product = new int[dim-1];

	    origin  = new RealVector ( dim );
	    tmp     = new RealVector ( dim );
	    tmp_    = new RealVector ( dim );

	    computeDelta();

	    VectorOperations.assign( c, -0.5 );
	}
    }

    void computeSizes() {
	
	computeRanges();

	for( int i=0; i<dim; i++ ) {
	    size[i] = range[i].max - range[i].min + 3;
	}
		
	volume = 1;
	for( int i=0; i<dim; i++ )
	    volume *= size[i];

	if( dim > 1 ) {
	    product[0] = size[0];
	    for( int i=1; i<dim-1; i++ )
		product[i] = product[i-1] * size[i];
	}

	for( int i=0; i<dim; i++ )
	    origin.re[i] = range[i].min - 1;

	originOffset = origin.re[0];

	for( int i=1; i<dim; i++ )
	    originOffset += product[i-1] * origin.re[i];
    }
    
    public final int index( RealVector y ) {
	return index( y.re );
    }

    public final int index( double [] yRe ) {

	double I = yRe[0];

	for( int i=1; i<dim; i++ )
	    I += product[i-1] * yRe[i];

	return (int)(I-originOffset+0.5);

    }

    void addLatticePoint( RealVector y ) {

	if( numOfLatticePoints + deltaLength >= maxNumOfLatticePoints){
	    maxNumOfLatticePoints *= 2 + deltaLength;
	    latticePoints.resize(maxNumOfLatticePoints, dim);
	}

	//debug
	//System.out.println( ( y.minus( new RealVector( c ) ).times( B.times(y.minus( new RealVector( c )  )) ) )  + " < " + sqrRadius );
	
	for (int i=0; i<deltaLength; i++ ) {
	    	    	    
	    tmp.assignPlus( y, delta[i] );
	  
	    final int indexOfYPlusOffset = index( tmp );

	    if( useSymmetrie ) {
		tmp_.assignNeg( tmp );
	    }
	    
	    if( !isSelected[ indexOfYPlusOffset ]  && 
		( !useSymmetrie || !isSelected[ index( tmp_ ) ] ) ) {
		
		latticePoints.setRow( numOfLatticePoints ++, tmp );
		isSelected[ indexOfYPlusOffset ] = true;
	    }      	      
	}	
    }


    public void update() {

	if( uptodate )
	    return;

	computeSizes();

	if( volume > isSelected.length )
	    isSelected = new boolean[ volume ];
	else
	    java.util.Arrays.fill( isSelected, false );

	tmp.assignZero();

	latticePoints.setRow( 0, tmp );

	isSelected[  index( tmp ) ] = true;

	numOfLatticePoints = 1;
  
	rho[dim-1] = sqrRadius;

	searchPoints( dim - 1 );

	uptodate = true;
    }

    double getMaxSqrRadiusOfDeltas() {
	double max = 0;

	for( int i=1; i<deltaLength; i++ ) {
	    tmp.assignTimes( B, delta[i] );
	    double aSqrRadius = tmp.dot( delta[i] );

	    if( aSqrRadius > max )
		max = aSqrRadius;
	}

	return max;
    }
}
