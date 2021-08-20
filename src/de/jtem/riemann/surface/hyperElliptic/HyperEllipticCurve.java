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

package de.jtem.riemann.surface.hyperElliptic;

import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexConstant;
import de.jtem.riemann.surface.AbstractAlgebraicCurveWithForms;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.surface.DistinguishedPoint;
import de.jtem.riemann.surface.Origin;
import de.jtem.riemann.surface.SingularPoint;

public class HyperEllipticCurve extends AbstractAlgebraicCurveWithForms {

    private static final long serialVersionUID = 1L;

    protected HyperEllipticCurve() {

    }

    public HyperEllipticCurve(BranchPoint bp[],
                              SingularPoint[] sp, DistinguishedPoint[] dp,
                              Complex origin) {

      branchPoint = bp;
      singularPoint = sp;
      distinguishedPoint = dp;

      this.origin = new Origin(origin.getRe(), origin.getIm());
    }

    public HyperEllipticCurve(BranchPoint bp[], Complex origin) {

      this( bp, null, null, origin );
    }

    public int  getNumOfSheets() {
	return 2;
    }

    private Complex tmpQ         = new Complex();
    private Complex tmpDQdLambda = new Complex();
    private Complex tmpDQdMu     = new Complex();

    private Complex tmpMu        = new Complex();
    private Complex tmpDif       = new Complex();
    private Complex tmpCoords    = new Complex();
    private Complex tmp          = new Complex();

    public void initCompute() {
	super.initCompute();

	tmpQ         = new Complex();
	tmpDQdLambda = new Complex();
	tmpDQdMu     = new Complex();

	tmpMu        = new Complex();
	tmpDif       = new Complex();
	tmpCoords    = new Complex();
	tmp          = new Complex();
    }

    private Complex lastMuAtOrigin;

    public void updateMuAtOrigin() {
	super.updateMuAtOrigin();

	if( lastMuAtOrigin == null )
	    lastMuAtOrigin = new Complex();

	if( lastMuAtOrigin.distSqr( musAtOrigin[1] ) <
	    lastMuAtOrigin.distSqr( musAtOrigin[0] ) ) {

	    Complex tmp    = musAtOrigin[0];

	    musAtOrigin[0] = musAtOrigin[1];
	    musAtOrigin[1] = tmp;
	}

	lastMuAtOrigin.assign( musAtOrigin[0] );

    }


    public void getMusAt( Complex lambda, ComplexVector mus ) {

	getProduct( lambda, tmpQ );

	tmpMu.assignSqrt( tmpQ );

	if( tmpMu.absSqr() < 1e-16 ) {
	    mus.newSize(1);

	    mus.set( 0, ComplexConstant.ZERO );
	} else {

	    mus.newSize( 2 );

	    mus.set( 0, tmpMu );

	    tmpMu.assignTimes(-1 );
	    mus.set( 1, tmpMu );
	}
    }


    public void getProduct( Complex lambda, Complex product ) {

	product.assign( 1 );

	for( int i=0; i<numOfBranchPoints; i++ ) {
		tmpDif.assignMinus( lambda, branchPoint[i].getCoords( tmpCoords ) );

		product.assignTimes( tmpDif );
	}
    }

    final void getProduct( Complex lambda, Complex product, int omittedBranchPoint ) {

	product.assign( 1, 0 );

	for( int i=0; i<numOfBranchPoints; i++ ) {

	    if( i != omittedBranchPoint ) {
		tmpDif.assignMinus( lambda, branchPoint[i].getCoords( tmpCoords ) );

		product.assignTimes( tmpDif );
	    }
	}
    }


    public void evaluateCurve( Complex lambda, Complex mu,
			       Complex Q, Complex dQdLambda, Complex dQdMu ) {

	getProduct( lambda, tmpQ );

	tmp.assignTimes( mu, mu );
	Q.assignMinus( tmpQ, tmp );

	dQdMu.assignTimes( mu, -2 );

	/* comptue dQdLambda */

	dQdLambda.assign( 0, 0 );

	if( tmpQ.absSqr() > 1e-6 ) {

	    for( int i=0; i<numOfBranchPoints; i++ ) {

		tmpDif.assignMinus( lambda, branchPoint[i].getCoords( tmpCoords ) );

		tmp.assignDivide( tmpQ, tmpDif );

		dQdLambda.assignPlus( tmp );
	    }
	} else {

	    for( int i=0; i<numOfBranchPoints; i++ ) {

		getProduct( lambda, tmp, i );

		dQdLambda.assignPlus( tmp );
	    }
	}
    }

    public int  getNumOfForms() {
	return  ( numOfBranchPoints - 1 ) / 2 + numOfSingularPoints;
    }

    public void getF( Complex lambda, Complex mu,
		      Complex Q, Complex dQdLambda, Complex dQdMu, ComplexVector f ) {

	tmp.assignInvert( mu );

	f.newSize( numOfForms );

	f.set( 0, tmp );

	for( int j = 1; j<numOfForms; j++ ) {

	    tmp.assignTimes( lambda );

	    f.set( j, tmp );
	}
    }

    public void getG( Complex lambda, Complex mu,
		      Complex Q, Complex dQdLambda, Complex dQdMu, ComplexVector g ) {

	tmp.assignInvert( dQdLambda );
	tmp.assignTimes( 2 );

	g.newSize( numOfForms );

	g.set( 0, tmp );

	for( int j = 1; j<numOfForms; j++ ) {

	    tmp.assignTimes( lambda );

	    g.set( j, tmp );
	}
    }

}

