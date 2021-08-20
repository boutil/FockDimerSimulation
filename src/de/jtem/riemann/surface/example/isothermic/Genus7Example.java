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

package de.jtem.riemann.surface.example.isothermic;

import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.polynomial.ComplexPolynomial;
import de.jtem.riemann.constrainedComplex.Dependent;
import de.jtem.riemann.constrainedComplex.Real;
import de.jtem.riemann.constrainedComplex.RealSymmetrie;
import de.jtem.riemann.surface.AbstractAlgebraicCurveWithForms;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.surface.Origin;
import de.jtem.riemann.surface.SurfacePoint;

public class Genus7Example extends AbstractAlgebraicCurveWithForms {
    
  private static final long serialVersionUID = 1L;

  final static double EPS4 = EPS*EPS*EPS*EPS;
  final static double EPS2 = EPS*EPS*EPS*EPS;

  ComplexPolynomial R = new ComplexPolynomial( 6 );
  ComplexPolynomial P = new ComplexPolynomial( 2 );
  ComplexPolynomial Q = new ComplexPolynomial( 6 );

  Complex[] rootsOfR = new Complex[6];
  Complex[] rootsOfQ;

  BranchPoint [] branchPointOfR;
  BranchPoint [] branchPointOfQ;

  public Genus7Example() {

    origin = new Origin( 0.1, 0.1 );

    P.assignMonomial( 2 );

    branchPoint    = new BranchPoint[ 12 ];
    branchPointOfR = new BranchPoint[  6 ];
    branchPointOfQ = new BranchPoint[  6 ];
      
    branchPoint[0] = branchPointOfR[0] = new BranchPoint( -2,  0 );
    branchPoint[1] = branchPointOfR[1] = new BranchPoint(  2,  0 );

    branchPoint[2] = branchPointOfR[2] = new BranchPoint(  0, -1 );
    branchPoint[3] = branchPointOfR[3] = new BranchPoint(  0,  1 );

    branchPoint[4] = branchPointOfR[4] = new BranchPoint(  1,  1 );
    branchPoint[5] = branchPointOfR[5] = new BranchPoint(  1, -1 );
   
    new Real( branchPoint[0] );
    new Real( branchPoint[1] );

    new RealSymmetrie( branchPoint[2], branchPoint[3] );
    new RealSymmetrie( branchPoint[4], branchPoint[5] );

    updatePolynomials();

    rootsOfQ = Q.getRoots();

    for( int i = 0; i < rootsOfQ.length; i++ ) {

      branchPoint[i+6] = branchPointOfQ[i] = new BranchPoint( rootsOfQ[i] );

      new Dependent( branchPointOfQ[i] ); 
    }

    //update();
    initCompute();
    update();
  }

  Complex lambda = new Complex();

  int [] permutation = new int[6];
  boolean [] isUsed = new boolean[6];

  void matchCoords( SurfacePoint [] point, Complex [] coords ) {

    for( int i=0; i<isUsed.length; i++ )
      isUsed[i] = false;

    for( int i=0; i<point.length; i++ ) {

      point[i].getCoords( lambda );

      double minDistSqr = Double.POSITIVE_INFINITY;

      for( int j=0; j<coords.length; j++ ) {

	if( !isUsed[j] ) {

	  final double distSqr = lambda.distSqr( coords[j] );

	  if( distSqr < minDistSqr ) {
	    permutation[i] = j;
	    minDistSqr = distSqr;
	  }
	}
      }

      isUsed[ permutation[i] ] = true;
    }

    for( int i=0; i<point.length; i++ )
      point[i].setCoords( coords[ permutation[ i ] ] );
  }

  void updateBranchPointsOfQ() {  // need to check continuity

    Q.getRoots( rootsOfQ );
    
    matchCoords( branchPointOfQ, rootsOfQ );
  }

  void updatePolynomials() {

    for( int i=0; i<branchPointOfR.length; i++ ) {
      rootsOfR[i] = branchPointOfR[i].getCoords( rootsOfR[i] );
    }

    R.setRoots( rootsOfR );

    Q.assignTimes( P, P );
    Q.assignMinus( R );
  }

  public void update() {

    if( uptodate )
      return;

    updatePolynomials();
    updateBranchPointsOfQ();

    super.update();


  }

  public int  getNumOfSheets() {
    return 4;
  }

  public void initCompute() {
    super.initCompute();
  }

  Complex pOfLambda = new Complex();
  Complex qOfLambda = new Complex();
  Complex rOfLambda = new Complex();
  Complex x         = new Complex();
  Complex y         = new Complex();
  Complex z         = new Complex();
  Complex muSqr     = new Complex();

  Complex pPrimeOfLambda = new Complex();
  Complex qPrimeOfLambda = new Complex();

  public void getMusAt( Complex lambda, ComplexVector mus ) {
    
    P.eval( lambda, pOfLambda );
    Q.eval( lambda, qOfLambda );
    
    rOfLambda.assignSqr( pOfLambda );
    rOfLambda.assignMinus( qOfLambda );

    final double absSqrOfROfLambda = rOfLambda.absSqr();
    final double absSqrOfQOfLambda = qOfLambda.absSqr();
    
    if( absSqrOfROfLambda < EPS4 ) {
      
      if( absSqrOfQOfLambda < EPS4 ) {
	mus.newSize( 1 );
	mus.set( 0, 0, 0 );
      } else {
	mus.newSize( 2 );
	x.assignSqrt( pOfLambda );  mus.set( 0, x );
	x.assignNeg();              mus.set( 1, x );
      }
    } else {

      if( absSqrOfQOfLambda < EPS2 ) {
	mus.newSize( 3 );
	mus.set( 0, 0, 0 );
	y.assignTimes( pOfLambda, 2 );
	x.assignSqrt( pOfLambda );  mus.set( 1, x );
	x.assignNeg();              mus.set( 2, x );
      } else {
 	mus.newSize( 4 );

	z.assignSqrt( rOfLambda );

	y.assignPlus( pOfLambda, z );
	x.assignSqrt( y );	    mus.set( 0, x );
	x.assignNeg(  x );          mus.set( 1, x );
	
	y.assignMinus( pOfLambda, z );
	x.assignSqrt( y );	    mus.set( 2, x );
	x.assignNeg(  x );          mus.set( 3, x );    
      }
    }


//     System.out.println( "check mus:" );

//     for( int i=0; i<mus.size(); i++ ) {
//       Field.Complex mu = mus.get(i);
//       Field.Complex s = new Field.Complex();

//       evaluateCurve( lambda, mu, s, new Field.Complex(), new Field.Complex() );

//       System.out.println( s );
//     }
    
  }
    
   

  public void evaluateCurve( Complex lambda, Complex mu, 
			     Complex s, Complex dQdLambda, Complex dQdMu ) {

    muSqr.assignSqr( mu );
 
    P.eval( lambda, pOfLambda );
    Q.eval( lambda, qOfLambda );
    
    s.assignTimes( pOfLambda, -2 );
    s.assignPlus(  muSqr );
    s.assignTimes( muSqr );
    s.assignPlus(  qOfLambda );

    P.evalDerivative( lambda, 1, pPrimeOfLambda );
    Q.evalDerivative( lambda, 1, qPrimeOfLambda );
    
    dQdLambda.assignTimes( pPrimeOfLambda, -2 );
    dQdLambda.assignTimes( muSqr );
    dQdLambda.assignPlus(  qPrimeOfLambda );

    dQdMu.assignMinus( muSqr, pOfLambda );
    dQdMu.assignTimes( mu );
    dQdMu.assignTimes( 4 );
  }


  public int  getNumOfForms() {
    return 7;
  }

  Complex w = new Complex();

  
  void computeForms( Complex lambda, Complex mu, 
		     Complex dSdLambda, Complex dSdMu, ComplexVector f ) {


    w.assign      ( mu );                f.set( 0, w );    
    w.assignTimes ( mu,     mu );        f.set( 1, w );
    w.assignTimes ( lambda, mu );        f.set( 2, w ); 


    w.assign      ( 1, 0 );              f.set( 3, w );
    w.assign      (         lambda );    f.set( 4, w );
    w.assignTimes ( w,      lambda );    f.set( 5, w );
    w.assignTimes ( w,      lambda );    f.set( 6, w );
  }

  public void getF( Complex lambda, Complex mu, 
		    Complex S, Complex dSdLambda, Complex dSdMu, ComplexVector f ) {
   
    computeForms( lambda, mu, dSdLambda, dSdMu, f );

    f.assignDivide( dSdMu );
  }    


  public void getG( Complex lambda, Complex mu, 
		    Complex S, Complex dSdLambda, Complex dSdMu, ComplexVector g ) {

    computeForms( lambda, mu, dSdLambda, dSdMu, g );

    g.assignDivide( dSdLambda );
    g.assignNeg();
  }



  public  BranchPoint [] getBranchPointsOfR() {
    return branchPointOfR;
  }

  public  BranchPoint [] getBranchPointsOfQ() {
    return branchPointOfQ;
  }

  public ComplexPolynomial getP() {
    return P;
  }

  public ComplexPolynomial getR() {
    return R;
  }

  public ComplexPolynomial getQ() {
    return Q;
  }

}

