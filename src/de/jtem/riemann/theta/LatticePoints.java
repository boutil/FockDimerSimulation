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

import de.jtem.blas.RealMatrix;
import de.jtem.blas.RealVector;
import de.jtem.numericalMethods.algebra.linear.decompose.Cholesky;

/**
 * 
 * @author schmies
 */
abstract class LatticePoints {

	protected final RealMatrix B;
	protected final RealMatrix [] T;
	protected final RealVector [] t;
	protected final RealVector [] abs_t;
	protected final RealVector N;
	protected final int dim;

	
	LatticePoints( int dim ) {
		this.dim = dim;
		
		B = new RealMatrix(dim);
		
		N = new RealVector(dim);
	
		T = new RealMatrix[dim+1];
		t = new RealVector[dim ];
		abs_t = new RealVector[dim ];
		
		C = new RealVector[dim+1];
		
		T[dim] = new RealMatrix(dim);
		C[dim] = new RealVector(dim);

		for( int g=1; g<dim; g++) {
			C[g] = new RealVector(g);
		}
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
	public void setB(RealMatrix  B) {
		
		if( B.getNumCols() != dim || B.getNumRows() != dim )
			throw new IllegalArgumentException( "Matrix has wrong dimension");
		
		uptodate = false;
				
		this.B.assign( B );
	
		Cholesky.decompose( B.re, T[dim].re );
	
		T[dim].assignTranspose();
		
		for( int g=1; g<dim; g++) {
			T[g] = T[dim].getBlock(0,0,g,g);
			
			RealVector v  = T[dim].getCol(g); v.setSize( g );
			t[g] = T[g].invert().times( v );
			
			abs_t[g] = new RealVector(g);
			for( int i=0; i<g; i++ ) {
				abs_t[g].re[i] = Math.abs(t[g].re[i] );
			}
		}
	}

	/**
	 * Returns center of ellipsoid.
	 */
	public final RealVector getCenter() {
		return new RealVector ( C[dim] );
	}

	/**
	 * Sets center of ellipsoid to <code>c</code>.
	 * @param c vector of length <code>dim</code>
	 */
	public final void setCenter(RealVector c) {
		if( c.size() != dim )
			throw new IllegalArgumentException( "center has wrong dimension" );
	
		C[dim].assign(c);
	
		uptodate = false;
	}

	protected static double sqr(double x) { return x*x; }

	public final double getRadius() {
		return r;
	}

	/** 
	 * Sets radius of ellipsoid.
	 * @param r radius of ellipsoid
	 */
	public final void setRadius(double  r) {
		this.r =r;
	
		uptodate = false;
	}

	protected final RealVector [] C;
	protected double r;
	protected boolean uptodate = false;

}
