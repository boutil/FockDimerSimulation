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

/** 
 * Iterates thru all lattice points inside an ellipsoids. 
 *
 * The ellipsoid is provided by a symmetric positive definite matrix 
 * <code>B</code>, a center vector <code>C</code> and, a radius <code>r</code>.
 *  Given this data this class determines all lattice points 
 * <code>N&isin;Z<sup>g</sup></code> that lay inside the specified ellipsoids, 
 * e.g. the lattice points fulfill:
 * <p align=center>
 *   <code>
 *    (N-c)<sup>t</sup> B (N-c) < r<sup>2</sup> with |c<sub>i</sub>-C<sub>i</sub>| < d<sub>i</sub>
 *   </code>,
 * </p>
 *
 * The iterator is not quite standard. It is optimized to determine the indices
 * for many ellipsoid with same the shape, e.g. the same matrix B, but for
 * different centers and radii. Here an example code:
 * 
 * <pre>
 *     double radius = 5;
 *		RealVector center = new RealVector(2);
 *		RealMatrix B = new RealMatrix( new double [][] {{2,1},{1,2}});
 *		
 *		LatticePointsInEllipsoidIterator lp 
 *			= new LatticePointsInEllipsoidIterator( B );
 *	
 *		lp.startIteration( radius, center);
 *			
 *		while(lp.hasNext() ) {
 *			System.out.println( new RealVector(lp.n) );
 *		}
 * </pre>		
 * @author Markus Schmies
 */

public class LatticePointsInEllipsoidIterator extends LatticePointsInEllipsoid implements Serializable, Cloneable {

	private static final long serialVersionUID = 1L;

	final int [] upper , lower;
	
	final double [] radius,  Tgg, n;
	final double [][] c;
	
	LatticePointsInEllipsoidIterator( int dim ) {
		super( dim );
		
		upper = new int[dim];
		lower = new int[dim];
		
		n = new double[dim];
		radius = new double[dim];
		Tgg = new double[dim+1];
		
		c = new double[dim+1][];
		for( int g=1; g<=dim; g++) {
			c[g] = new double[ g ];
		}
		
	}
	
	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in zero and has radius 1.
	 * @param B symmetric and positive definite.
	 */
	public LatticePointsInEllipsoidIterator( RealMatrix B ) {
		this( B, 1 );
	}

	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in zero and has radius <code>r</code>.
	 * @param B symmetric and positive definite.
	 * @param r radius of ellipsoid
	 */
	public LatticePointsInEllipsoidIterator( RealMatrix B, double r ) {	
		this( B, new RealVector(B.getNumCols()), r );
	}

	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in <code>c</code> and has radius <code>r</code>.
	 * @param B symmetric and positive definite.
	 * @param C center of ellipsoid
	 * @param r radius of ellipsoid
	 */
	public LatticePointsInEllipsoidIterator( RealMatrix B, RealVector C, double r ) {	
		this( B.getNumCols() );
		setB( B );
		setCenter( C );
		setRadius( r );
	}
	
	/**
	 * Sets matrix specifiying the ellipsoid.
	 * @param B symmetric and positive definite
	 */
	public final void setB(RealMatrix  B ) {
		super.setB(B);
	
		Tgg[dim] = T[dim].get(dim-1,dim-1);
		
		for( int g=1; g<dim; g++) {
			Tgg[g]=T[g].get(g-1,g-1);
		}
	}
	
	/**
	 * Iterates to next lattice point, which is returned in N.
	 * @param N next lattice point on output
	 * @return true if it has a next lattice point
	 */
	public boolean hasNext( RealVector  N) {
		boolean out = hasNext();
		N.assign(n);
		return out;
	}
	
	boolean hasNext() {
		n[0]++;
		if( n[0] <= upper[0]  )
			return true;
		if( dim == 1)
			return false;
		
		return hasNext(2);
	}
	
	boolean hasNext( int g ) {
		final int h = g-1;
		
		n[h]++;

		if( n[h] <= upper[h] ) {
			final double [] t_h = t[h].re;
			
			for( int j=0; j<h; j++  ) {
					c[h][j] -= t_h[j];
			}
				
			radius[h-1] = Math.sqrt( sqr( radius[h] ) - sqr( Tgg[g] * (n[h]-c[g][h]) ) );
			startIteration( h  );
			
			return true;
		}
		
		if( g==dim )
			return false;
		
		return hasNext( g+1 );
	}
	
	/**
	 * Starts iteration thru all lattice points enclosed by the ellipsoid
	 * prescribed by the shape induced by the matrix B, the
	 * radius r, and center C.
	 */
	public void startIteration() {	
		startIteration(r, C[dim].re);
	}
	
	/**
	 * Starts iteration thru all lattice points enclosed by the ellipsoid
	 * prescribed by the shape induced by the matrix B and the
	 * radius r, center C.
	 * @param r radius
	 * @param C center
	 */
	public void startIteration( double r, RealVector C  ) {	
		startIteration(r, C.re);
	}
	
	void startIteration( double r, double [] c ) {
		System.arraycopy( c, 0, this.c[dim], 0, c.length );
		radius[dim-1] = r;
		startIteration( dim );
		n[0]--;
	}
	
	void startIteration( int g ) {
		
		final int h = g-1;
		
		final double Tgg =  this.Tgg[g];
		
		final double r = radius[h];
		
		lower[h] = (int)Math.floor( c[g][h]  - r / Tgg + 1 );
		upper[h] = (int)Math.floor( c[g][h] + r / Tgg  );

		double n_g = n[h]=lower[h];
				
		if( h==0 )
			return;
		
		final double c_g = c[g][h];
		
		final double [] t_h = t[h].re;
		
		for( int j=0; j<h; j++  ) {
			c[h][j] = c[g][j] - t_h[j]*( n_g - c_g );
		}
				
		radius[h-1] = Math.sqrt( r*r - sqr( Tgg * (n_g-c_g) ) );
		
		startIteration( h  );
	}
	
	
	public static void main( String [] arg ) {
		double radius = 5;
		RealVector center = new RealVector(2);
		RealMatrix B = new RealMatrix( new double [][] {{2,1},{1,2}});
		
		LatticePointsInEllipsoidIterator lp 
			= new LatticePointsInEllipsoidIterator( 
					B, center, radius );
	
		lp.startIteration( radius, center);
			
		while(lp.hasNext() ) {
			System.out.println( new RealVector(lp.n) );
		}
		
		System.out.println( new LatticePointsInEllipsoid( new RealMatrix( new double [][] {{2,1},{1,2}}), center, radius  ).getLatticePoints() );
	}
	
}
