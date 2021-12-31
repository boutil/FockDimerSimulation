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

import java.io.Serializable;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.blas.IntegerMatrix;
import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.calculus.integration.NewtonCotes;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.surface.DistinguishedPoint;
import de.jtem.riemann.surface.Origin;
import de.jtem.riemann.surface.RamifiedCovering;
import de.jtem.riemann.surface.SingularPoint;
import de.jtem.riemann.surface.homologie.CanonicalBasis;
import de.jtem.riemann.surface.homologie.Transform;

public class HyperEllipticSurface extends RamifiedCovering implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    protected int genus;
    protected int numOfForms;

    protected Complex lastMu;
    protected Complex thisMu;
    protected Complex lastLambda;
    protected Complex thisLambda;
    protected Complex muAtOrigin;

    ComplexVector integral;
    ComplexMatrix discretizer;

    int numOfSteps;

    Complex tmp, tmpDif;

    boolean muAtOriginUptodate = false;

    ComplexVector [] branchPointIntegral;
    ComplexVector [] singularPointIntegral;
    ComplexVector [] distinguishedPointIntegral;

    boolean surfacePointIntegralsUptodate = false;

    ComplexVector [] branchPointLoop;
    ComplexVector [] singularPointLoop;
    ComplexVector [] distinguishedPointPath;

    ComplexMatrix branchPointIntegralMatrix;
    ComplexMatrix singularPointIntegralMatrix;
    ComplexMatrix distinguishedPointIntegralMatrix;

    ComplexVector edgeIntegralVector;

    boolean symmetrizePeriodMatrix = true;

    protected CanonicalBasis cb;

    public HyperEllipticSurface( BranchPoint bp [],
				 SingularPoint [] sp, DistinguishedPoint [] dp, Complex origin ) {

	branchPoint = bp;
	singularPoint = sp;
	distinguishedPoint = dp;

	origin = new Origin( origin.re, origin.im );

	update();
    }

    public Transform getTransform() {
	return getHomologieTransform();
    }

    public Transform getHomologieTransform() {
	return cb;
    }

    public void initCompute() {

	super.initCompute();

	cb = new CanonicalBasis( getMonodromyMatrix(),
				 getSingularMatrix(), getDistinguishedVector() );

	genus      = ( numOfBranchPoints - 1 ) / 2;

	numOfForms = genus + numOfSingularPoints;

	branchPointIntegral       = new ComplexVector[numOfBranchPoints];
	branchPointIntegralMatrix = new ComplexMatrix( numOfBranchPoints, 2 );


	for( int i=0; i<numOfBranchPoints; i++ )
	    branchPointIntegral[i] = new ComplexVector( numOfForms );

	branchPointLoop = new ComplexVector[numOfBranchPoints ];

	if( numOfSingularPoints > 0 ) {

	    singularPointIntegral       = new ComplexVector[numOfSingularPoints];
	    singularPointIntegralMatrix = new ComplexMatrix( numOfSingularPoints, 2 );

	    for( int i=0; i<numOfSingularPoints; i++ )
		singularPointIntegral[i] = new ComplexVector( numOfForms );

	    singularPointLoop = new ComplexVector[ numOfSingularPoints ];

	} else {
	    singularPointIntegral = null;
	    singularPointLoop = null;
	}

	if( numOfDistinguishedPoints > 0 ) {

	    distinguishedPointIntegral = new ComplexVector[numOfDistinguishedPoints];

	    for( int i=0; i<numOfDistinguishedPoints; i++ )
		distinguishedPointIntegral[i] = new ComplexVector( numOfForms );

	    distinguishedPointPath           = new ComplexVector[ numOfDistinguishedPoints ];
	    distinguishedPointIntegralMatrix = new ComplexMatrix( numOfDistinguishedPoints, 2 );


	} else {
	    distinguishedPointIntegral = null;
	    distinguishedPointPath = null;
	}

	edgeIntegralVector = new ComplexVector( cb.getNumOfEdges() );

	discretizer = new ComplexMatrix( numOfSteps, numOfForms );
	integral    = new ComplexVector( numOfForms );

	/* we need to fix a sheet */
	lastMu      = new Complex();
	thisMu      = new Complex();

    	lastLambda  = new Complex();
    	thisLambda  = new Complex();

	tmpDif      = new Complex();
	tmp         = new Complex();

	muAtOrigin  = new Complex();

	muAtOriginUptodate = false;
    }

    public final boolean getSurfacePointIntegralsUptodate() {
	return surfacePointIntegralsUptodate;
    }

    public Complex getMuAtOrigin() {
	return new Complex( muAtOrigin );
    }

    public void updateMuAtOrigin() {

	if( muAtOriginUptodate )
	    return;

	thisMu.assign( muAtOrigin );

	muAtOrigin.assign( computeMu( origin ) );

	muAtOriginUptodate = true;
    }


    public void outdate() {
	super.outdate();

	 muAtOriginUptodate = surfacePointIntegralsUptodate = false;
    }

    public void update() {

	if( uptodate )
	    return;

	super.update();

	updateMuAtOrigin();

	computeLoops();
    }

    /**
       * Get the value of genus.
       * @return Value of genus.
       */
    public final int getGenus() {return genus;}

    /**
       * Get the value of numOfForms.
       * @return Value of numOfForms.
       */
    public final int getNumOfForms() {return numOfForms;}

    public final boolean getSymmetrizePeriodMatrix() {
	return symmetrizePeriodMatrix;
    }

    public void setSymmetrizePeriodMatrix( boolean aBool ) {
	if( aBool == symmetrizePeriodMatrix )
	    return;

	symmetrizePeriodMatrix = aBool;
	outdate();
    }
    /**
       * Get the value of integral.
       * @return Value of integral.
       */
    public ComplexVector getIntegral() {return integral;}

    /**
       * Set the value of integral.
       * @param v  Value to assign to integral.
       */
    public void setIntegral(ComplexVector  v) {this.integral = v;}

    /**
       * Get the value of discretizer.
       * @return Value of discretizer.
       */
    public ComplexMatrix getDiscretizer() {return discretizer;}

    /**
       * Set the value of discretizer.
       * @param v  Value to assign to discretizer.
       */
    public void setDiscretizer(ComplexMatrix  v) {
	this.discretizer = v;

	checkDiscretizerCapacity();
    }

    /**
       * Get the value of numOfSteps.
       * @return Value of numOfSteps.
       */
    public int getNumOfSteps() {return numOfSteps;}

    /**
       * Set the value of numOfSteps.
       * @param v  Value to assign to numOfSteps.
       */
    public void setNumOfSteps(int  v) {
	this.numOfSteps = v;
	outdate();
    }


    public ComplexVector[] getBranchPointLoop() {
	return branchPointLoop;
    }

    public ComplexVector[] getSingularPointLoop() {
	return singularPointLoop;
    }

    public ComplexVector[] getDistinguishedPointPath() {
	return distinguishedPointPath;
    }

    public IntegerMatrix getMonodromyMatrix() {
	IntegerMatrix monodrom = new IntegerMatrix ( numOfBranchPoints, 2 );

	for( int i=0; i<numOfBranchPoints; i++ ) {
	    monodrom.set( i, 0, 1 );
	    monodrom.set( i, 1, 0 );
	}

	return monodrom;
    }

    public IntegerMatrix getSingularMatrix() {

	if( numOfSingularPoints < 1 )
	    return null;

	IntegerMatrix aMatrix= new IntegerMatrix ( numOfSingularPoints, 2 );

	for( int i=0; i<numOfSingularPoints; i++ ) {
	    aMatrix.set( i, 0, 1 );
	    aMatrix.set( i, 1, 1 );
	}

	return aMatrix;
    }

    public final void hyperEllipticSort( Complex lastValue, Complex thisValue ) {
	if( lastValue.re * thisValue.re + lastValue.im * thisValue.im < 0 ) {
	    thisValue.re *= -1;
	    thisValue.im *= -1;
	}
    }

    public final int hyperEllipticSortFactor( final Complex lastValue, final Complex thisValue ) {
	return lastValue.re * thisValue.re + lastValue.im * thisValue.im < 0 ? -1 : 1;
    }

    public Complex computeMu( Complex at ) {


	computeMu( at, -1 );

	hyperEllipticSort( lastMu, thisMu );

	return thisMu;
    }

    public Complex computeMu( final Complex at, final int excludedBranchPointNumber ) {

	lastLambda.assign( thisLambda );
	thisLambda.assign( at );

	lastMu.assign( thisMu );
	thisMu.assign( 1, 0 );

	for( int i=0; i<numOfBranchPoints; i++ )
	    if( i != excludedBranchPointNumber ) {

		tmpDif.assignMinus( at, branchPoint[i].getCoords( tmp ) );

		thisMu.assignTimes( tmpDif );
	    }

	thisMu.assignSqrt();

	return thisMu;
    }


    private void checkDiscretizerCapacity() {
	checkDiscretizerCapacity( numOfForms );
    }

    private void checkDiscretizerCapacity( int numForms ) {

	if( discretizer.getNumRows() < numForms || discretizer.getNumCols() != numOfSteps )
	    discretizer.newSize( numForms, numOfSteps );
    }

    private Complex value  = new Complex();

    public void evaluateFormsAt( final Complex lambda,
				 final Complex mu, final ComplexVector values ) {

	final int numForms = values.size();

	value.assignInvert( mu );

	values.set( 0, value );

	for( int j = 1; j<numForms; j++ ) {

	    value.assignTimes( lambda );

	    values.set( j, value );
	}
    }

    public void evaluateFormsAt( final Complex lambda, final ComplexVector values ) {

	computeMu( lambda );

	evaluateFormsAt( lambda, thisMu, values );
    }

    public void evaluateFormsAtBranchPoint( final Complex Q, final ComplexVector values ) {

	int indexOfQ = getIndexOfBranchPointWithCoords( Q, 1e-20 );

	if( indexOfQ == - 1 )
	    throw new IllegalArgumentException
		("distance of " + Q + " to any branch point is larger than 1e-10!" );

	computeMu( Q, indexOfQ );

	evaluateFormsAt( Q, thisMu, values );
    }

    private ComplexVector values = new ComplexVector(0);

    private Complex V         = new Complex();
    private Complex lambda    = new Complex();
    private Complex QMinusP   = new Complex();
    private Complex lastValue = new Complex();

    private void discretizeForms( Complex P, Complex Q, int numForms ) {

	if( numForms == 0 )
	    return;

	values.newSize( numForms );

	V.assignMinus( Q, P );
	V.assignDivide( numOfSteps-1 );

	lambda.assign( P );

	checkDiscretizerCapacity();

	updateMuAtOrigin();

	for( int i = 0; i<numOfSteps; i++, lambda.assignPlus( V ) ) {

	    evaluateFormsAt( lambda, values );

	    discretizer.setCol( i, values );
	}
    }

    private void discretizeFormsIntoBranchPoint( Complex P, Complex Q, int numForms ) {

	if( numForms == 0 )
	    return;

	checkDiscretizerCapacity();

	updateMuAtOrigin();

	QMinusP.assignMinus( Q, P );

	for( int i = 0; i<numOfSteps - 1; i++ ) {

	    double t = Math.PI / 2 * (double)i / (numOfSteps-1);

	    lambda.assignTimes( QMinusP, Math.sin( t ) );
	    lambda.assignPlus( P );

	    evaluateFormsAt( lambda, values );
	    values.assignTimes( Math.cos(t) );
	    discretizer.setCol( i, values );
	}

	values.get( 0, lastValue );

	evaluateFormsAtBranchPoint( Q, values );

	tmp.assignDivide( QMinusP, -2 );
	tmp.assignSqrt();

	values.assignDivide( tmp );

	values.get( 0, value );

	values.assignTimes( hyperEllipticSortFactor( lastValue, value ) );

	discretizer.setCol( numOfSteps - 1, values );

	//discretizer.getCol( numOfSteps - 2 ).print("1");
	//discretizer.getCol( numOfSteps - 1 ).print("2 ");
    }

  private void discretizeFormsIntoInfinity( Complex P, int numForms ) {

	if( numForms == 0 )
	    return;

	checkDiscretizerCapacity();

	updateMuAtOrigin();

	for( int i = 0; i<numOfSteps - 1; i++ ) {

	    double t = 1 - i / (numOfSteps-1.);

	    lambda.assignTimes( P, 1/t );

	    computeMu( lambda );

	    value.assignInvert( thisMu );
	    value.assignTimes( 1/t/t );

	    discretizer.set( 0, i, value );

	    for( int j = 1; j<numForms; j++ ) {
		value.assignTimes( lambda );
		discretizer.set( j, i, value );
	    }
	}

	for( int j=0; j<genus-1; j++ )
	    discretizer.set( j, numOfSteps-1, 0, 0 );

	value.assignInvert( P );
	value.assignTimes( value );

	hyperEllipticSort( discretizer.get( genus-1, numOfSteps-2 ), value );

	discretizer.set( genus-1, numOfSteps-1, value );

	for( int j = genus; j<numForms; j++ )
	    discretizer.set( j, numOfSteps-1, Double.MAX_VALUE, Double.MAX_VALUE );
    }

    public final void integrateForms( Complex P, Complex Q, int numForms ) {
	discretizeForms( P, Q, numForms );
	tmp.assignMinus( Q, P );
	sumDiscretizer( numForms, tmp );
    }

    public final void integrateFormsIntoBranchPoint( Complex P, Complex Q, int numForms ) {
	discretizeFormsIntoBranchPoint( P, Q, numForms );
	tmp.assignMinus( Q, P );
	tmp.assignTimes( Math.PI / 2 );
	sumDiscretizer( numForms, tmp );
    }

    public final void integrateFormsIntoInfinity( Complex P, int numForms ) {
	discretizeFormsIntoInfinity( P, numForms );
	sumDiscretizer( numForms, P );
    }

    private final void sumDiscretizer( int numForms, Complex V ) {

	Complex value = new Complex();

	if( numForms == 0 )
	    return;

	if( integral.size() < numForms )
	    integral.newSize( numForms );

	for( int j=0; j<numForms; j++ ) {

	    value.assign( NewtonCotes.highestOrderSum( discretizer.re[j] ),
		       NewtonCotes.highestOrderSum( discretizer.im[j] ) );

	    value.assignTimes( V );

	    integral.set( j, value );
	}
    }

    public final void integrateForms( Complex P, Complex Q, ComplexVector aVec, int numForms ) {

	ComplexVector saveIntegral = integral;

	aVec.newSize( numForms );
	setIntegral( aVec );

	integrateForms( P, Q, numForms );

	integral = saveIntegral;
    }

    public final void integrateFormsIntoBranchPoint( Complex P, Complex Q,
						     ComplexVector aVec, int numForms ) {

	ComplexVector saveIntegral = integral;

	aVec.newSize( numForms );
	setIntegral( aVec );

	integrateFormsIntoBranchPoint( P, Q, numForms );

	integral = saveIntegral;
    }

    public final void integrateFormsIntoInfinity( Complex P, ComplexVector aVec, int numForms ) {

	ComplexVector saveIntegral = integral;

	aVec.newSize( numForms );
	setIntegral( aVec );

	integrateFormsIntoInfinity( P, numForms );

	integral = saveIntegral;
    }


    public void integrateForms( Complex P, Complex Q, ComplexVector aVec ) {
	integrateForms( P, Q, aVec, numOfForms );
    }

    public void integrateFormsIntoBranchPoint( Complex P, Complex Q, ComplexVector aVec ) {
	integrateFormsIntoBranchPoint( P, Q, aVec, numOfForms );
    }

    public void integrateFormsIntoInfinity( Complex P, ComplexVector aVec ) {
	integrateFormsIntoInfinity( P, aVec, numOfForms );
    }


    public void integrateHolomorphicForms( Complex P, Complex Q, ComplexVector aVec ) {
	integrateForms( P, Q, aVec, genus );
     }

    public final void integrateHolomorphicFormsIntoBranchPoint( Complex P, Complex Q, ComplexVector aVec ) {
	integrateFormsIntoBranchPoint( P, Q, aVec, genus );
    }

    public final void integrateHolomorphicFormsIntoInfinity( Complex P, ComplexVector aVec ) {
	integrateFormsIntoInfinity( P, aVec, genus );
    }


    public final void integrateForms( ComplexVector path, ComplexVector aVec, int numForms ) {
	integrateForms( path, aVec, numForms, false );
    }

    public final void integrateFormsIntoBranchPoint( ComplexVector path, ComplexVector aVec, int numForms ) {
	integrateForms( path, aVec, numForms, true );
    }

    private final Complex P = new Complex();
    private final Complex Q = new Complex();

    public final void integrateForms( ComplexVector path, ComplexVector aVec,
				int numForms, boolean lastIsBranchPoint ) {

	int numOfPoints = path.size();

	aVec.newSize( numForms );
	aVec.assignZero();

	if( numOfPoints < 2 )
          return;

	for( int i = 1; i<numOfPoints; i++ ) {

	    path.get( i-1, P );
	    path.get( i,   Q );

	    if( i == numOfPoints-1 && lastIsBranchPoint )
		integrateFormsIntoBranchPoint( P, Q, numForms );
	    else
		integrateForms( P, Q, numForms );

	    for( int j = 0; j<numForms; j++ ) {
		aVec.re[j] += integral.re[j];
		aVec.im[j] += integral.im[j];
	    }
	}
    }

    public final void integrateForms( ComplexVector path, ComplexVector aVec ) {
	integrateForms( path, aVec, numOfForms, false );
    }

    public final void integrateFormsIntoBranchPoint( ComplexVector path, ComplexVector aVec ) {
	integrateForms( path, aVec, numOfForms, true );
    }

    public final void integrateHolomorphicForms( ComplexVector path, ComplexVector aVec ) {
	integrateForms( path, aVec, genus, false );
    }

    public final void integrateHolomorphicFormsIntoBranchPoint( ComplexVector path, ComplexVector aVec ) {
	integrateForms( path, aVec, genus, true );
    }

    public void computeLoops() {

	computeLoops( branchPointLoop, singularPointLoop, distinguishedPointPath );
    }

    void computeSurfacePointIntegrals( int numForms ) {

	update();

	if( surfacePointIntegralsUptodate )
	    return;

	surfacePointIntegralsUptodate = true;

	for( int i=0; i<numOfBranchPoints; i++ ) {
	    thisMu.assign( muAtOrigin );
	    integrateForms( branchPointLoop[i], branchPointIntegral[i], numForms );
	}

	if( numOfSingularPoints > 0 ) {

	    for( int i=0; i<numOfSingularPoints; i++ ) {
		thisMu.assign( muAtOrigin );
		integrateForms( singularPointLoop[i], singularPointIntegral[i], numForms );
	    }
	}

	if( numOfDistinguishedPoints > 0 ) {

	    for( int i=0; i<numOfDistinguishedPoints; i++ ) {
		thisMu.assign( muAtOrigin );
		integrateForms( distinguishedPointPath[i],
				distinguishedPointIntegral[i], numForms,
				distinguishedPoint[i].isAlsoBranchPoint  );
	    }
	}
    }


    public ComplexMatrix getEdgeIntegrals( int numForms ) {

	ComplexMatrix edgeIntegrals = new ComplexMatrix( cb.getNumOfEdges(), genus );

	computeEdgeIntegrals( edgeIntegrals, numForms );

	return edgeIntegrals;
    }

    public void computeEdgeIntegrals( ComplexMatrix edgeIntegrals, int numForms ) {
	computeSurfacePointIntegrals( numForms );

	edgeIntegrals.newSize( cb.getNumOfPaths(), numForms );

	for( int i=0; i<numForms; i++ ) {

	    for( int j=0; j<numOfBranchPoints; j++ ) {
		double rr = branchPointIntegral[j].re[i];
		double ii = branchPointIntegral[j].im[i];

		branchPointIntegralMatrix.re[j][0] =  rr;
		branchPointIntegralMatrix.im[j][0] =  ii;
		branchPointIntegralMatrix.re[j][1] = -rr;
		branchPointIntegralMatrix.im[j][1] = -ii;
	    }

	    if( numOfSingularPoints > 0 ) {

		for( int j=0; j<numOfSingularPoints; j++ ) {
		    double rr = singularPointIntegral[j].re[i];
		    double ii = singularPointIntegral[j].im[i];

		    singularPointIntegralMatrix.re[j][0] =  rr;
		    singularPointIntegralMatrix.im[j][0] =  ii;
		    singularPointIntegralMatrix.re[j][1] = -rr;
		    singularPointIntegralMatrix.im[j][1] = -ii;
		    }
	    }

	    if( numOfDistinguishedPoints > 0 ) {

		for( int j=0; j<numOfDistinguishedPoints; j++ ) {
		    double rr = distinguishedPointIntegral[j].re[i];
		    double ii = distinguishedPointIntegral[j].im[i];

		    distinguishedPointIntegralMatrix.re[j][0] =  rr;
		    distinguishedPointIntegralMatrix.im[j][0] =  ii;
		    distinguishedPointIntegralMatrix.re[j][1] = -rr;
		    distinguishedPointIntegralMatrix.im[j][1] = -ii;
		}
	    }

	    cb.computeEdgeIntegrals( branchPointIntegralMatrix,
						     singularPointIntegralMatrix,
						     distinguishedPointIntegralMatrix, edgeIntegralVector );

	    edgeIntegrals.setCol( i, edgeIntegralVector );
	}

    }

    public ComplexMatrix getAPeriods() {

	ComplexMatrix edgeIntegrals = getEdgeIntegrals( genus );

	return cb.getAPeriodsForEdgeIntegrals( edgeIntegrals );


    }

    public ComplexMatrix getBPeriods() {

	ComplexMatrix edgeIntegrals = getEdgeIntegrals( genus );

	return cb.getBPeriodsForEdgeIntegrals( edgeIntegrals );
    }

    public static void symmetrizePeriodMatrix( ComplexMatrix periodMatrix ) {

	int genus = periodMatrix.getNumRows();

	for( int i=0; i<genus; i++ )
	    for( int j=i+1; j<genus; j++ ) {
		periodMatrix.re[i][j] = ( periodMatrix.re[i][j] + periodMatrix.re[j][i] )/2;
		periodMatrix.im[j][i] = ( periodMatrix.im[i][j] + periodMatrix.im[j][i] )/2;

		periodMatrix.re[j][i] =   periodMatrix.re[i][j];
		periodMatrix.im[j][i] =   periodMatrix.im[i][j];

	    }
    }

    public ComplexMatrix getPeriodMatrix() {

	ComplexMatrix periodMatrix = new ComplexMatrix( genus );

	computePeriodMatrix( periodMatrix );

	return periodMatrix;
    }

    public void computePeriodMatrix( ComplexMatrix periodMatrix ) {

	ComplexMatrix edgeIntegrals = getEdgeIntegrals( genus );

	cb.computePeriodMatrixForEdgeIntegrals( edgeIntegrals, periodMatrix );

	if( symmetrizePeriodMatrix )
	    symmetrizePeriodMatrix( periodMatrix );
    }
}




