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
 * Determines all lattice points inside the union of all ellipsoids of a certain shape
 * and radius but with different centers, which all lie in a prescribed g-dimensional interval.
 *
 * The ellipsoid is provided by a symmetric positive definite matrix 
 * <code>B</code>, a center vector <code>C</code> and, a radius <code>r</code> and
 * a deviation  vector <code>D</code>.  Given this data this class determines all lattice points 
 * <code>N&isin;Z<sup>g</sup></code> that lay inside the specified ellipsoids, 
 * e.g. the lattice points fulfill:
 * <p align=center>
 *   <code>
 *    (N-c)<sup>t</sup> B (N-c) < r<sup>2</sup> with |c<sub>i</sub>-C<sub>i</sub>| < d<sub>i</sub>
 *   </code>,
 * </p>
 *
 * The lattice points are provided as rows of a matrix (see: 
 * {@link #getLatticePoints()}).
 *
 * @author Markus Schmies
 */

public class LatticePointsForUniformApproximation extends LatticePointsList implements Serializable, Cloneable {

	private static final long serialVersionUID = 1L;

	final RealVector [] D;
	
	LatticePointsForUniformApproximation( int dim ) {
		super( dim );
		
		D = new RealVector[dim+1];
		
		D[dim] = new RealVector(dim);
		
		for( int g=1; g<dim; g++) {
			D[g] = new RealVector(g);
		}
	}
	
	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in zero and has radius 1.
	 * @param B symmetric and positive definite.
	 */
	public LatticePointsForUniformApproximation( RealMatrix B ) {
		this( B, 1 );
	}

	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in zero and has radius <code>r</code>.
	 * @param B symmetric and positive definite.
	 * @param r radius of ellipsoid
	 */
	public LatticePointsForUniformApproximation( RealMatrix B, double r ) {	
		this( B, new RealVector(B.getNumCols()), 0.5, r );
	}

	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in <code>c</code> and has radius <code>r</code>.
	 * @param B symmetric and positive definite.
	 * @param C center of ellipsoid
	 * @param r radius of ellipsoid
	 */
	public LatticePointsForUniformApproximation( RealMatrix B, RealVector C, double r ) {	
		this( B, C, 0.5, r );
	}

	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in <code>c</code> and has radius <code>r</code>.
	 * @param B symmetric and positive definite.
	 * @param C center of ellipsoid
	 * @param D deviation of center of ellipsoid
	 * @param r radius of ellipsoid
	 */
	public LatticePointsForUniformApproximation( RealMatrix B, RealVector C, double D, double r ) {	
		this( B.getNumCols() );
		setB( B );
		setCenter( C );
		setRadius( r );
		setDeviation( D );
	}
	
	/**
	 * Creates object that delivers the lattice points inside ellipsoid
	 * given by the matrix <code>B</code>.
	 * The ellipsoid is centered in <code>c</code> and has radius <code>r</code>.
	 * @param B symmetric and positive definite.
	 * @param C center of ellipsoid
	 * @param D deviation of center of ellipsoid
	 * @param r radius of ellipsoid
	 */
	public LatticePointsForUniformApproximation( RealMatrix B, RealVector C, RealVector D, double r ) {	
		this( B.getNumCols() );
		setB( B );
		setCenter( C );
		setRadius( r );
		setDeviation( D );
	}
	
	/**
	 * Sets deviation range of center of ellipsoid to <code>c</code>.
	 * @param d vector of length <code>dim</code> with non-negative entries.
	 */
	public final void setDeviation( double  d ) {
	
		D[dim].assign(d);

		uptodate = false;
	}
	

	/**
	 * Sets deviation range of center of ellipsoid to <code>c</code>.
	 * @param c vector of length <code>dim</code>
	 */
	public final void setDeviation( RealVector d ) {
		if( d.size() != dim )
			throw new IllegalArgumentException( "deviation has wrong dimension" );

		D[dim].assign(d);

		uptodate = false;
	}
	
	private void searchPoints( int g, double r, RealVector n) {
	
		final int h = g-1;
		
		final double Tgg = T[g].get(h,h);
		
		final int lower = (int)Math.floor( C[g].re[h] - D[g].re[h] - r / Tgg + 1 );
		final int upper = (int)Math.floor( C[g].re[h] + D[g].re[h] + r / Tgg  );

		if( g==1 ) {
//			if( maxDeviation < D[g].re[h] )
//				maxDeviation =D[g].re[h];
//			
			for( int i=lower; i<=upper; i++ ) {
				n.re[h]=i;
				this.addLatticePoint(n);
			}
		} else {

			final double c_g = C[g].re[h];
			final double d_g = D[g].re[h];
		
			final double [] t_h = t[h].re;
			final double [] abs_t_h = abs_t[h].re;
			
			for( int j=0; j<h; j++  ) {
				D[h].re[j] = D[g].re[j] + abs_t_h[j]*d_g;
			}
			
			for( int i=lower; i<=upper; i++ ) {
				double n_g = n.re[g-1]=i;
				
				for( int j=0; j<h; j++  ) {
					C[h].re[j] = C[g].re[j] - t_h[j]*( n_g - c_g );
				}
				
				// minSqr is min ( n_g - x )^2 with cg + d_g <= x <= c_g + d_g
				double minSqr =  
					n_g-c_g<-d_g ? sqr(n_g-c_g+d_g) : 
						n_g-c_g> d_g ? sqr(n_g-c_g-d_g) : 0;
						
						searchPoints( g-1,  Math.sqrt( r*r - Tgg*Tgg * minSqr ), n );
			}		
		}
	}
	
	protected void update() {
		
		if( uptodate )
			return;
		
		numOfLatticePoints = 0;
		//maxDeviation = 0;
		
		searchPoints( dim, r, N );
		
		uptodate = true;
		
		//System.out.println( "maxDeviation =" + maxDeviation );
	}
	
	public static void main( String [] arg ) {		
		System.out.println( new LatticePointsForUniformApproximation( new RealMatrix( new double [][] {{7.5,0},{0,191}}), 5.8 ).getLatticePoints() );
	}
	
}




