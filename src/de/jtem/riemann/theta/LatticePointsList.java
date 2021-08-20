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

abstract class LatticePointsList extends  LatticePoints implements Serializable, Cloneable {

	protected final RealMatrix latticePoints;
	protected int numOfLatticePoints = 0;
	int maxNumOfLatticePoints = 100;
	
	LatticePointsList( int dim ) {
		super(dim);
		
		latticePoints = new RealMatrix ( 100, dim );
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

	abstract void update();
	
	/**
	 * Returns matrix containing the lattice points as string.
	 * @see #getLatticePoints()
	 */
	public final String toString() {
		return getLatticePoints().toString();
	}

	protected void addLatticePoint(RealVector y) {
		
		if(numOfLatticePoints == maxNumOfLatticePoints){
			maxNumOfLatticePoints *= 2;
			latticePoints.resize(maxNumOfLatticePoints, dim);
		}
		
		System.arraycopy( y.re, 0, latticePoints.re[numOfLatticePoints ++], 0, dim );
		
		//latticePoints.setRow(numOfLatticePoints ++, y );
	}

}
	