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
import de.jtem.numericalMethods.algebra.linear.decompose.Cholesky;

/** 
 * Determines all lattice points inside a specified ellipsoid.
 *
 * The ellipsoid is provided by a symmetric positive definite matrix 
 * <code>B</code>, a center vector <code>C</code> and, a radius <code>r</code>.
 * Given this data this class determines all lattice points 
 * <code>N&isin;Z<sup>g</sup></code> that lay inside the specified ellipsoid, 
 * e.g. the lattice points fulfill:
 * <p align=center>
 *   <code>
 *    (N-C)<sup>t</sup> B (N-C) < r<sup>2</sup> 
 *   </code>,
 * </p>
 *
 * The lattice points are provided as rows of a matrix (see: 
 * {@link #getLatticePoints()}).
 *
 * @author Markus Schmies
 */

public class LatticePointsInEllipsoid_old implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    boolean uptodate = false;

    RealMatrix B = new RealMatrix();

    RealVector  y;

    double [] rho;

    double [][] v;
    double [][] l;

    double [] c;     

    RealMatrix latticePoints;

    int	numOfLatticePoints = 0;

    int maxNumOfLatticePoints = 100;

    int dim;

    double sqrRadius;
    
    Range [] range;   // Range is inner class


    /**
     * Creates object that delivers the lattice points inside ellipsoid
     * given by the matrix <code>B</code>.
     * The ellipsoid is centered in zero and has radius 1.
     * @param B symmetric and positive definite.
     */
    public LatticePointsInEllipsoid_old( RealMatrix B ) {
	this( B, 1 );
    }

    /**
     * Creates object that delivers the lattice points inside ellipsoid
     * given by the matrix <code>B</code>.
     * The ellipsoid is centered in zero and has radius <code>r</code>.
     * @param B symmetric and positive definite.
     * @param r radius of ellipsoid
     */
    public LatticePointsInEllipsoid_old( RealMatrix B, double r ) {	

	dim = -1;
	setB( B );
	setRadius( r );
    }

    /**
     * Creates object that delivers the lattice points inside ellipsoid
     * given by the matrix <code>B</code>.
     * The ellipsoid is centered in <code>c</code> and has radius <code>r</code>.
     * @param B symmetric and positive definite.
     * @param C center of ellipsoid
     * @param r radius of ellipsoid
     */
    public LatticePointsInEllipsoid_old( RealMatrix B, RealVector C, double r ) {	

	dim = -1;
	setB( B );
	setCenter( C );
	setRadius( r );
    }

    /**
     * Returns dimension of vectorspace.
     */
    public final int getDim() {return dim;}
    
    void setDim(int  newDim) { 

	if( newDim != dim ) {
	    dim = newDim;

	    l = new double[dim][dim];

	    latticePoints = new RealMatrix ( 100, dim );

	    y = new RealVector ( dim );

	    c   = new double[ dim ];	    
	    rho = new double[ dim ];
	
	    v = new double[dim][];
	    
	    range = new Range[ dim ];

	    for( int i=0; i<dim; i++ ) {
		range[i] = new Range(0,0);

		v[i] = new double[ i + 1 ];		
	    }
	}
    }

    /**
     * Returns lattice points in ellipsoid.
     * The lattice points are the rows of the returned matrix.
     * @return matrix containing the lattice points
     */
    public final RealMatrix getLatticePoints() {
	update();

	RealMatrix croppedField = new RealMatrix ( latticePoints );
	
	croppedField.resize( numOfLatticePoints, dim );

	return croppedField;
    }

    /**
     * Returns number of lattice points in specified ellipsoid.
     */
    public final int getNumOfLatticePoints() {
	update();
	return numOfLatticePoints;
    }

    /**
     * Returns matrix specifying the ellipsoid.
     */
    public final RealMatrix getB() {
	return new RealMatrix (B);
    }
    
    /**
     * Sets matrix specifiying the ellipsoid.
     * @return B symmetric and positive definite
     */
    public final void setB(RealMatrix  B ) {
	
	uptodate = false;
	
	this.B.assign( B );

	setDim( B.getNumRows() );

	Cholesky.decompose( B.re, l );

    }
    
    /**
     * Returns center of ellipsoid.
     */
    public final RealVector getCenter() {
	return new RealVector ( c );
    }
    
    /**
     * Sets center of ellipsoid to <code>c</code>.
     * @param c vector of length <code>dim</code>
     */
    public final void setCenter( RealVector c ) {
	if( c.size() != dim )
	    throw new IllegalArgumentException( "center has wrong dimension" );


	VectorOperations.assign( this.c, c.re );

	uptodate = false;
    }

    final void setSqrRadius(double  newSqrRadius) {

	uptodate = false;

	sqrRadius = newSqrRadius;
    }

    void update() {

	if( uptodate )
	    return;

	numOfLatticePoints = 0;

	rho[dim-1] = sqrRadius;

	searchPoints( dim-1 );

	uptodate = true;
    }
 
    /** 
     * Returns radius of ellipsoid.
     */
    public final double getRadius() {
	return Math.sqrt( sqrRadius );
    }
    
    /** 
     * Sets radius of ellipsoid.
     * @param r radius of ellipsoid
     */
    public final void setRadius( double  r ) {
	setSqrRadius( r * r );
    }

    /**
     * Returns matrix containing the lattice points as string.
     * @see #getLatticePoints()
     */
    public final String toString() {
	return getLatticePoints().toString();
    }

    static final class Range implements java.io.Serializable {

	private static final long serialVersionUID = 1L;

        int min, max;

	Range( double aMin, double aMax ) {
	    set( aMin, aMax );
	}

	Range set( double aMin, double aMax ) {
	    min = (int)Math.floor( aMin + 1 ); 
	    max = (int)Math.floor( aMax     ); 

	    return this;
	} 

	Range set( int aMin, int aMax ) {
	    min = aMin;
	    max = aMax;

	    return this;
	} 

	Range set( Range otherRange ) {
	    min = otherRange.min;
	    max = otherRange.max;

	    return this;
	} 
    }

    void addLatticePoint( RealVector y ) {
	
	if(numOfLatticePoints == maxNumOfLatticePoints){
	    maxNumOfLatticePoints *= 2;
	    latticePoints.resize(maxNumOfLatticePoints, dim);
	}
	
	latticePoints.setRow(numOfLatticePoints ++, y );
    }

    void searchPoints( int  k ) {

	final double vkk = v[k][k];
	final double lkk = l[k][k];

	final double root  = Math.sqrt( rho[k] / lkk / lkk );

	final double major = -vkk / lkk + c[k];

	range[k].set( major - root, major + root );

	//System.out.println( "range " + k + " [ " + range[k].min + " , " + range[k].max + " ] " );

	final double min = range[k].min;
	final double max = range[k].max + 0.5;

	for( double x = min; x < max; x += 1 ) {
	    y.re[k] = x;
	    
	    if( k>0 ) {
		final double xMinusCk = x - c[k];
		final double tmp = xMinusCk * lkk + vkk;

		rho[k-1] = rho[k] - tmp * tmp;
				
		for( int i=0; i<k; i++ )
		    v[k-1][i] = v[k][i] + xMinusCk * l[k][i];

		searchPoints( k-1 );
	    } else {
		addLatticePoint ( y );
	    }
	}
    }    

    final void computeRanges() {

	RealMatrix B_inv = B.invert();
	
	for (int i=0; i<dim; i++ ) {
	    
	    double max = Math.sqrt( sqrRadius * B_inv.re[i][i] );
	    
	    range[i].set( -max + c[i], max + c[i] );
	}	
    }

    final boolean isSelectedLatticePoint( RealVector x ) {
	return isSelectedLatticePoint ( x, latticePoints, numOfLatticePoints);
    }

    static boolean isSelectedLatticePoint( RealVector x, RealMatrix field, int fieldLength ) {

	for( int i=0; i<fieldLength; i++ )
	    if( x.equals( field.getRow(i) ) ) {
		return true;
	    }
	return false;
    }
    
    private int checkLength, checkVolume;

    void check(){

	checkLength = checkVolume = 0;

	computeRanges();

	check( 0 );

	System.out.println("fieldLength = "+ numOfLatticePoints+
			   "  checkLength = " + checkLength+
			   "  checkVolume = " + checkVolume );    
    }

    private void check( int k ) {
	if( k != dim )  {
	    for (int i = range[k].min; i <= range[k].max; i++ ) {
		y.re[k] = i;
		check(k+1);
	    }
	}
	
	else {
	    
	    double product = B.times( y.minus( new RealVector(c) ) ).dot( y.minus( new RealVector(c) ) );

	    if( product <= sqrRadius  ) {

		checkVolume++;

		if( !isSelectedLatticePoint ( y ) ) {

		    checkLength++;

		    System.out.println("Is not in field (" + product + ") " + y );
		}
	    }
	}
    }

}




