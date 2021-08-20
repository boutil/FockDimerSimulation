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

package de.jtem.riemann.surface;

import java.io.Serializable;

import de.jtem.blas.ComplexVector;
import de.jtem.blas.IntegerMatrix;
import de.jtem.blas.IntegerVector;
import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.calculus.odeSolving.Extrap;
import de.jtem.numericalMethods.calculus.odeSolving.ODE;
import de.jtem.numericalMethods.calculus.odeSolving.ODEIntermediateResultListener;
import de.jtem.riemann.analyticContinuation.AlgebraicCurve;
import de.jtem.riemann.surface.homologie.CanonicalBasis;
import de.jtem.riemann.surface.homologie.Transform;

public abstract class AbstractAlgebraicCurve extends RamifiedCoveringVariation implements AlgebraicCurve, Serializable, Cloneable {

	private static final long serialVersionUID = 1L;

	protected int genus;
	//      protected int numOfForms;
	protected int numOfSheets;

	final Complex tmp = new Complex();
	final Complex tmpDif = new Complex();

	protected Complex[] musAtOrigin;
	protected Complex[] musAtDistinguishedPoints;

	//      ComplexVector tmpIntegral, tmpVector, tmpEdgeIntegralVector;

	ComplexVector musVector;

	//      ComplexMatrix [] branchPointIntegral;
	//      ComplexMatrix [] singularPointIntegral;
	//      ComplexMatrix [] distinguishedPointIntegral;

	ComplexVector[] branchPointLoop;
	ComplexVector[] singularPointLoop;
	ComplexVector[] distinguishedPointPath;

	//      ComplexMatrix branchPointIntegralMatrix;
	//      ComplexMatrix singularPointIntegralMatrix;
	//      ComplexMatrix distinguishedPointIntegralMatrix;

	//      ComplexVector edgeIntegralVector;

	//      protected ComplexMatrix edgeIntegrals;

	//      boolean symmetrizePeriodMatrix = true;

	//      FormFunction ff;

	//      Extrap integrator;

	//      double [] I, H;

	//    double eps = 1e-8;

	//      protected ComplexMatrix aPeriods;
	//      protected ComplexMatrix bPeriods;
	//      protected ComplexMatrix formTransform;
	//      protected ComplexMatrix periodMatrix;

	protected CanonicalBasis cb;

	IntegerMatrix monodromyMatrix;
	IntegerMatrix singularMatrix;
	IntegerVector distinguishedVector;

	public abstract int getNumOfSheets();
	public abstract void getMusAt(Complex lambda, ComplexVector mus);

	/* (non-Javadoc)
	 * @see de.jtem.riemann.surface.AlgebraicCurve#evaluateCurve(de.jtem.mfc.field.Complex, de.jtem.mfc.field.Complex, de.jtem.mfc.field.Complex, de.jtem.mfc.field.Complex, de.jtem.mfc.field.Complex)
	 */
	public abstract void evaluateCurve(Complex lambda, Complex mu, Complex Q, Complex dQdLambda, Complex dQdMu);

	//      public abstract int  getNumOfForms();
	//      public abstract void getF( Field.Complex lambda, Field.Complex mu,
	//  			       Field.Complex Q, Field.Complex dQdLambda, Field.Complex dQdMu, ComplexVector f );
	//      public abstract void getG( Field.Complex lambda, Field.Complex mu,
	//  			       Field.Complex Q, Field.Complex dQdLambda, Field.Complex dQdMu, ComplexVector g );

	public Transform getTransform() {
		return getHomologieTransform();
	}

	public Transform getHomologieTransform() {
		return cb;
	}

	public void initHomologie() {
		cb = new CanonicalBasis(getMonodromyMatrix(), getSingularMatrix(), getDistinguishedVector());
	}

	public void initCompute() {

		super.initCompute();

		//  numOfForms  = getNumOfForms();
		numOfSheets = getNumOfSheets();

		//  	branchPointIntegral       = new ComplexMatrix[numOfBranchPoints];
		//  	branchPointIntegralMatrix = new ComplexMatrix( numOfBranchPoints, numOfSheets );

		//  	for( int i=0; i<numOfBranchPoints; i++ )
		//  	    branchPointIntegral[i] = new ComplexMatrix( numOfSheets, numOfForms );

		branchPointLoop = new ComplexVector[numOfBranchPoints];

		if (numOfSingularPoints > 0) {

			//  	    singularPointIntegral       = new ComplexMatrix[numOfSingularPoints];
			//  	    singularPointIntegralMatrix = new ComplexMatrix( numOfSingularPoints, numOfSheets );

			//  	    for( int i=0; i<numOfSingularPoints; i++ )
			//  		singularPointIntegral[i] = new ComplexMatrix( numOfSheets, numOfForms );

			singularPointLoop = new ComplexVector[numOfSingularPoints];

		} else {
			//      singularPointIntegral = null;
			singularPointLoop = null;
		}

		musAtDistinguishedPoints = new Complex[numOfDistinguishedPoints];

		if (numOfDistinguishedPoints > 0) {

			//  	    distinguishedPointIntegral = new ComplexMatrix[numOfDistinguishedPoints];

			for (int i = 0; i < numOfDistinguishedPoints; i++) {
				//  		distinguishedPointIntegral[i] = new ComplexMatrix( numOfForms, numOfForms );

				musAtDistinguishedPoints[i] = new Complex();
			}

			distinguishedPointPath = new ComplexVector[numOfDistinguishedPoints];
			//  	    distinguishedPointIntegralMatrix = new ComplexMatrix( numOfDistinguishedPoints,
			//  								  numOfSheets );

		} else {
			//  	    distinguishedPointIntegral = null;
			distinguishedPointPath = null;
		}

		musAtOrigin = new Complex[numOfSheets];

		musVector = new ComplexVector(numOfSheets);
		//  	tmpIntegral = new ComplexVector( numOfForms );
		//  	tmpVector   = new ComplexVector( numOfForms );

		/* stuff for the integration */
		//  	ff         = new FormFunction();

		//  	integrator = new Extrap( 4 + numOfForms * 2 );

		//  	I          = new double[ 5 + numOfForms * 2 ];
		//  	H          = new double[ 1 ];

		initHomologie();

		//  	tmpEdgeIntegralVector = new ComplexVector( cb.getNumOfPaths() );

		genus = cb.getGenus();

		//  	aPeriods           = new ComplexMatrix( genus );
		//  	bPeriods           = new ComplexMatrix( genus );
		//  	formTransform      = new ComplexMatrix( genus );
		//  	periodMatrix       = new ComplexMatrix( genus );

		//  	edgeIntegralVector = new ComplexVector( cb.getNumOfEdges() );
		//  	edgeIntegrals      = new ComplexMatrix( cb.getNumOfEdges(), numOfForms );
	}

	public Complex[] getMusAtOrigin() {
		return musAtOrigin;
	}

	public void updateMuAtOrigin() {

		getMusAt(origin, musVector);

		musVector.toArray(musAtOrigin);

		/// !!! need to match the new mus with the old
	}

	public void update() {

		if (uptodate)
			return;

		super.update();

		updateMuAtOrigin();
		updateLoops();
		//  	updateSurfacePointIntegrals();
		//  	updateEdgeIntegrals();
		//  	updateMatrices();
	}

	//      public double getEps() {
	//  	return eps;
	//      }

	//      public void setEps(double  v) {
	//  	eps = v;
	//  	outdate();
	//  	update();
	//      }

	public final int getGenus() {
		return genus;
	}

	//      public final boolean getSymmetrizePeriodMatrix() {
	//  	return symmetrizePeriodMatrix;
	//      }

	//      public void setSymmetrizePeriodMatrix( boolean aBool ) {
	//  	if( aBool == symmetrizePeriodMatrix )
	//  	    return;

	//  	symmetrizePeriodMatrix = aBool;
	//  	outdate();
	//  	update();
	//      }

	public ComplexVector[] getBranchPointLoop() {
		return branchPointLoop;
	}

	public ComplexVector[] getSingularPointLoop() {
		return singularPointLoop;
	}

	public ComplexVector[] getDistinguishedPointPath() {
		return distinguishedPointPath;
	}

	private class MonodromyFunction implements ODE {

		double t = 0;

		final Complex lambda = new Complex();
		final Complex mu = new Complex();

		final Complex lambdaDirection = new Complex();

		final Complex Q = new Complex();
		final Complex dQdLambda = new Complex();
		final Complex dQdMu = new Complex();

		final Complex dLambdaDt = new Complex();
		final Complex dMuDt = new Complex();

		final MonodromyCorrector mcf = new MonodromyCorrector();

		int numOfEvaluations = 0;

		private class MonodromyCorrector implements ODEIntermediateResultListener {

			double t;

			final Complex lambda = new Complex();
			final Complex mu = new Complex();
			

			final Corrector cf = new Corrector(AbstractAlgebraicCurve.this);

			public void intermediateResult(double t, double[] x, int offsetX) {
				//t = x[offsetX];

				lambda.assign(x[0 + offsetX], x[1 + offsetX]);
				mu.assign(x[2 + offsetX], x[3 + offsetX]);

				cf.correctMu(lambda, mu);

				x[0 + offsetX] = lambda.re;
				x[1 + offsetX] = lambda.im;

				x[2 + offsetX] = mu.re;
				x[3 + offsetX] = mu.im;

			}

		}

		public int getNumberOfEquations() {
			return 4;
		}

		public void eval(double t, double[] x, double[] y) {

			lambda.assign(x[0], x[1]);
			mu.assign(x[2], x[3]);

			evaluateCurve(lambda, mu, Q, dQdLambda, dQdMu);

			dLambdaDt.assign(lambdaDirection);

			dMuDt.assignDivide(dQdLambda, dQdMu);
			dMuDt.assignTimes(-1);
			dMuDt.assignTimes(lambdaDirection);

			y[0] = lambdaDirection.re;
			y[1] = lambdaDirection.im;

			y[2] = dMuDt.re;
			y[3] = dMuDt.im;

			numOfEvaluations++;
		}

	}

	public IntegerMatrix getMonodromyMatrix(ComplexVector[] loops) {

		MonodromyFunction mf = new MonodromyFunction();

		IntegerMatrix monodrom = new IntegerMatrix(loops.length, numOfSheets);

		Extrap odeSolver = new Extrap(4);

		updateMuAtOrigin();
		updateLoops();

		double[] Y = new double[5];
		double[] H = new double[1];

		for (int j = 0; j < loops.length; j++)
			for (int k = 0; k < numOfSheets; k++) {

				Y[1] = origin.re;
				Y[2] = origin.im;
				Y[3] = musAtOrigin[k].re;
				Y[4] = musAtOrigin[k].im;

				H[0] = 0.1;

				for (int i = 1; i < loops[j].size(); i++) {

					Y[0] = 0;

					mf.lambdaDirection.assignMinus(loops[j].get(i), loops[j].get(i - 1));

					odeSolver.odex(mf, Y, 1, H, 1e-7, 1e-7, mf.mcf);
				}

				Complex finalMu = new Complex(Y[3], Y[4]);

				int i;
				for (i = 0; i < numOfSheets; i++)
					if (finalMu.dist(musAtOrigin[i]) < 1e-8 * finalMu.abs())
						break;

				if (i == numOfSheets) {

					System.out.println("computed mu = " + finalMu);
					for (i = 0; i < numOfSheets; i++)
						System.out.println("possible match = " + musAtOrigin[i]);
					throw new RuntimeException("could not match mu in origin: ");
				}
				monodrom.set(j, k, i);
			}

		// System.out.println( "number evaluations: " + mf.numOfEvaluations );

		return monodrom;
	}

	public IntegerMatrix getMonodromyMatrix() {

		if (monodromyMatrix == null)
			monodromyMatrix = getMonodromyMatrix(branchPointLoop);

		return monodromyMatrix;
	}

	public IntegerMatrix getSingularMatrix() {

		if (singularMatrix == null && numOfSingularPoints != 0)
			singularMatrix = getMonodromyMatrix(singularPointLoop);

		return singularMatrix;
	}

	//      private class FormFunction implements DoubleArrayValued,
	//  					  DoubleArrayParametrized, Serializable, Cloneable {

	//  	private static final long serialVersionUID = 1L;

	//  	boolean integrateFormsIntoBranchPoint = false;

	//  	double  t = 0;

	//  	final Field.Complex lambda          = new Field.Complex();
	//  	final Field.Complex mu              = new Field.Complex();

	//  	final Field.Complex direction       = new Field.Complex();
	//  	final Field.Complex last            = new Field.Complex();

	//  	final Field.Complex Q               = new Field.Complex();
	//  	final Field.Complex dQdLambda       = new Field.Complex();
	//  	final Field.Complex dQdMu           = new Field.Complex();

	//  	final Field.Complex dLambdaDt       = new Field.Complex();
	//  	final Field.Complex dMuDt           = new Field.Complex();

	//  	final FormCorrector fcf       = new FormCorrector();

	//  	final ComplexVector values    = new ComplexVector( numOfForms );

	//  	int numOfEvaluations          = 0;

	//  	double localRadiusSqrOfLast   = 0;

	//  	private class FormCorrector implements DoubleArrayValued,
	//  					       DoubleArrayParametrized, Serializable, Cloneable {

	//  	    private static final long serialVersionUID = 1L;

	//  	    double t;

	//  	    final Field.Complex lambda      = new Field.Complex();
	//  	    final Field.Complex mu          = new Field.Complex();

	//  	    final Corrector cf        = new Corrector();

	//  	    public int getDoubleArrayParameterLength() {
	//  		return 5;
	//  	    }

	//  	    public void setByParameter( double [] p ) {
	//  		setDoubleArrayParameter ( p, 0 );
	//  	    }

	//  	    public void setDoubleArrayParameter( double [] p, int offset ) {
	//  		t = p[offset];

	//  		lambda.set( p[1+offset], p[2+offset] );
	//  		mu.    set( p[3+offset], p[4+offset] );

	//  		if( integrateFormsIntoBranchPoint )
	//  		    cf.correctLambda( lambda, mu );
	//  		else
	//  		    cf.correctMu( lambda, mu );
	//  	    }

	//  	    public double [] getDoubleArrayValue() {
	//  		double[] value = new double[5];

	//  		getDoubleArrayValue ( value, 0 );

	//  		return value;
	//  	    }

	//  	    public int getDoubleArrayValueLength() {
	//  		return 5;
	//  	    }

	//  	    public void getValue( double [] value ) {
	//  		getDoubleArrayValue ( value, 0 );
	//  	    }

	//  	    public void getDoubleArrayValue( double [] value, int offset ) {

	//  		value[0+offset] = t;

	//  		value[1+offset] = lambda.re;
	//  		value[2+offset] = lambda.im;

	//  		value[3+offset] = mu.re;
	//  		value[4+offset] = mu.im;
	//  	    }
	//  	}

	//  	public int getDoubleArrayParameterLength() {
	//  	    return 5 + numOfForms * 2;
	//  	}

	//  	public void setByParameter( double [] p ) {
	//  	    setDoubleArrayParameter ( p, 0 );
	//  	}

	//  	public void setDoubleArrayParameter( double [] p, int offset ) {
	//  	    t = p[offset];
	//  	    lambda.set( p[1+offset], p[2+offset] );
	//  	    mu.    set( p[3+offset], p[4+offset] );

	//  	    if( integrateFormsIntoBranchPoint ) {

	//  		fcf.cf.correctLambda( lambda, mu );

	//  		if( false && last.distSqr( lambda ) > localRadiusSqrOfLast ) {
	//  		    System.out.println( "localRadiusSqrOfQ = " + localRadiusSqrOfLast );
	//  		    System.out.println( "dist = " + last.distSqr( lambda ) );
	//  		    System.out.println( "lambda = " + lambda );
	//  		    System.out.println( "branch point = " + last );
	//  		    System.out.println( "t = " + t );

	//  		    // this is a big problem
	//  		    throw new RuntimeException( "exit local region" );
	//  		}

	//  		evaluateCurve( lambda, mu, Q, dQdLambda, dQdMu );

	//  		dMuDt.assign( direction );

	//  		dLambdaDt.assignDivide( dQdMu, dQdLambda );
	//  		dLambdaDt.assignTimes( -1 );
	//  		dLambdaDt.assignTimes( dMuDt );

	//  		getG( lambda, mu, Q, dQdLambda, dQdMu, values );
	//  	    } else {

	//  		fcf.cf.correctMu( lambda, mu );

	//  		evaluateCurve( lambda, mu, Q, dQdLambda, dQdMu );

	//  		dLambdaDt.assign( direction );

	//  		dMuDt.assignDivide( dQdLambda, dQdMu );
	//  		dMuDt.assignTimes( -1 );
	//  		dMuDt.assignTimes( dLambdaDt );

	//  		getF( lambda, mu, Q, dQdLambda, dQdMu, values );
	//  	    }

	//  	    values.assignTimes( direction );

	//  	    numOfEvaluations++;
	//  	}

	//  	public double [] getDoubleArrayValue() {
	//  	    double[] value = new double[5+numOfForms * 2];

	//  	    getDoubleArrayValue ( value, 0 );

	//  	    return value;
	//  	}

	//  	public int getDoubleArrayValueLength() {
	//  	    return 5 + numOfForms * 2;
	//  	}

	//  	public void getValue( double [] value ) {
	//  	    getDoubleArrayValue ( value, 0 );
	//  	}

	//  	public void getDoubleArrayValue( double [] value, int offset ) {

	//  	    if( integrateFormsIntoBranchPoint ) {
	//  		value[1+offset] = dLambdaDt.re;
	//  		value[2+offset] = dLambdaDt.im;

	//  		value[3+offset] = direction.re;
	//  		value[4+offset] = direction.im;
	//  	    } else {
	//  		value[1+offset] = direction.re;
	//  		value[2+offset] = direction.im;

	//  		value[3+offset] = dMuDt.re;
	//  		value[4+offset] = dMuDt.im;
	//  	    }

	//  	    for( int i=0, j=5; i<numOfForms; i++ ) {
	//  		value[j++ + offset] = values.re[i];
	//  		value[j++ + offset] = values.im[i];
	//  	    }
	//  	}
	//      }

	//      public final void integrateForms( ComplexVector path,
	//  				      Field.Complex muAtStartPoint,
	//  				      Field.Complex muAtEndPoint, ComplexVector aVec ) {
	//  	integrateForms( path, muAtStartPoint, muAtEndPoint, aVec, false );
	//      }

	//      public final void integrateFormsIntoBranchPoint( ComplexVector path,
	//  						     Field.Complex muAtStartPoint,
	//  						     Field.Complex muAtEndPoint,
	//  						     ComplexVector aVec ) {
	//  	integrateForms( path, muAtStartPoint, muAtEndPoint, aVec, true );
	//      }

	//      protected final void startIntegration( Field.Complex lambda, Field.Complex mu ) {

	//  	I[0] = 0;
	//  	I[1] = lambda.re;
	//  	I[2] = lambda.im;
	//  	I[3] = mu.re;
	//  	I[4] = mu.im;

	//  	for( int i=0, j=5; i<numOfForms; i++ ) {
	//  		 I[j++] = 0.0;
	//  		 I[j++] = 0.0;
	//  	}

	//  	H[0] = 0.1;

	//      }

	//      private final Field.Complex P = new Field.Complex();
	//      private final Field.Complex Q = new Field.Complex();

	//      protected final void integrateForms( Field.Complex Q ) {
	//  	ff.integrateFormsIntoBranchPoint = false;

	//  	P.set( I[1], I[2] );

	//  	double distQP = Q.dist(P);

	//  	if( distQP < 1e-10 )
	//  	    return;

	//  	ff.direction.assignMinus( Q, P );
	//  	ff.direction.assignDivide( distQP );

	//  	integrator.odex( ff, ff, I, I[0] + distQP, H, eps, eps, ff.fcf );
	//      }

	//      protected final void integrateForms( Field.Complex Q, ComplexVector aVec ) {

	//  	integrateForms( Q );

	//  	abelValue( aVec );
	//      }

	//      protected final void integrateForms( Field.Complex Q, Field.Complex muAtQ, ComplexVector aVec ) {

	//  	integrateForms( Q );

	//  	muAtQ.set( I[3], I[4] );

	//  	abelValue( aVec );
	//      }

	//      private final void abelValue( ComplexVector aVec ) {

	//  	aVec.newSize( numOfForms );

	//  	for( int i=0, j=5; i<numOfForms; i++ ) {
	//  	    aVec.set(i, I[j++], I[j++] );
	//  	}
	//      }

	//      public final void integrateForms( ComplexVector path,
	//  				      Field.Complex muAtStartPoint,
	//  				      Field.Complex muAtEndPoint,
	//  				      ComplexVector aVec, boolean lastIsBranchPoint ) {

	//  	int numOfPoints = path.size();

	//  	if( numOfPoints < 2 )
	//            return;

	//  	path.get( 0, P );

	//  	startIntegration( P, muAtStartPoint );

	//  	for( int i = 1; i<numOfPoints - 1; i++ ) {
	//  	    path.get(i,   Q );

	//  	    integrateForms( Q );
	//  	}

	//  	if( lastIsBranchPoint ) {

	//  	    path.get( numOfPoints-1, ff.last );

	//  	    BranchPoint lastPoint = getBranchPointWithCoords( ff.last );

	//  	    ff.localRadiusSqrOfLast  = lastPoint.getRadius() * 1.2;
	//  	    ff.localRadiusSqrOfLast *= ff.localRadiusSqrOfLast;

	//  	    ff.integrateFormsIntoBranchPoint = true;
	//  	    ff.direction.assign( -I[3], -I[4] );

	//  	    path.print( "path" );

	//  	    H[0] = 0.1;
	//  	    I[0] = 0;

	//  	    integrator.odex( ff, ff, I, 1, H, eps, eps, ff.fcf );

	//  	} else {
	//  	    path.get( numOfPoints-1, Q );

	//  	    integrateForms( Q );
	//  	}

	//  	muAtEndPoint.set( I[3], I[4] );

	//  	abelValue( aVec );

	//  	//	if( lastIsBranchPoint )
	//  	//    aVec.print( "to Z" );

	//      }

	public void updateLoops() {

		computeLoops(branchPointLoop, singularPointLoop, distinguishedPointPath);
	}

	//      void updateSurfacePointIntegrals() {

	//  	for( int i=0; i<numOfBranchPoints; i++ )
	//  	    for( int j=0; j<numOfSheets; j++ ) {

	//  		integrateForms( branchPointLoop[i], musAtOrigin[j], tmp, tmpVector );

	//  		branchPointIntegral[i].setRow( j, tmpVector );
	//  	    }

	//  	if( numOfSingularPoints > 0 )
	//  	    for( int i=0; i<numOfSingularPoints; i++ )
	//  		for( int j=0; j<numOfSheets; j++ ) {

	//  		    integrateForms( singularPointLoop[i], musAtOrigin[j], tmp, tmpVector );

	//  		    singularPointIntegral[i].setRow( j, tmpVector );
	//  		}

	//  	if( numOfDistinguishedPoints > 0 )

	//  	    for( int i=0; i<numOfDistinguishedPoints; i++ )
	//  		for( int j=0; j<numOfSheets; j++ ) {

	//  		    integrateForms( distinguishedPointPath[i],
	//  				    musAtOrigin[j], tmp, tmpVector,
	//  				    distinguishedPoint[i].isAlsoBranchPoint );

	//  		    if( j == distinguishedPoint[i].sheetNumber )
	//  			musAtDistinguishedPoints[i].set( tmp );

	//  		    distinguishedPointIntegral[i].setRow( j, tmpVector );
	//  		}
	//      }

	//      public ComplexMatrix getEdgeIntegrals() {
	//  	return edgeIntegrals;
	//      }

	//      public void updateEdgeIntegrals() {

	//  	edgeIntegrals.newSize( cb.getNumOfPaths(), numOfForms );

	//  	for( int i=0; i<numOfForms; i++ ) {

	//  	    for( int j=0; j<numOfBranchPoints; j++ )
	//  		for( int k=0; k<numOfSheets; k++ ) {

	//  		    branchPointIntegral[j].get( k, i, tmp );

	//  		    branchPointIntegralMatrix.set( j, k, tmp );
	//  		}

	//  	    if( numOfSingularPoints > 0 ) {

	//  		for( int j=0; j<numOfSingularPoints; j++ )
	//  		    for( int k=0; k<numOfSheets; k++ ) {

	//  		    singularPointIntegral[j].get( k, i, tmp );

	//  		    singularPointIntegralMatrix.set( j, k, tmp );
	//  		}
	//  	    }

	//  	    if( numOfDistinguishedPoints > 0 ) {
	//  		for( int j=0; j<numOfDistinguishedPoints; j++ )
	//  		    for( int k=0; k<numOfSheets; k++ ) {

	//  		    distinguishedPointIntegral[j].get( k, i, tmp );

	//  		    distinguishedPointIntegralMatrix.set( j, k, tmp );
	//  		}

	//  	    }

	//  	    cb.computeEdgeIntegrals( branchPointIntegralMatrix,
	//  				     singularPointIntegralMatrix,
	//  				     distinguishedPointIntegralMatrix,
	//  				     edgeIntegralVector );

	//  	    edgeIntegrals.setCol( i, edgeIntegralVector );
	//  	}
	//      }

	//      public void computeAPeriodsOfForm( int indexOfForm, ComplexVector edgeIntegralsOfForm ) {

	//  	edgeIntegrals.getCol( indexOfForm, tmpEdgeIntegralVector );

	//  	cb.computeAPeriodsForEdgeIntegrals( edgeIntegralsOfForm, tmpEdgeIntegralVector );
	//      }

	//      public void computeBPeriodsOfForm( int indexOfForm, ComplexVector edgeIntegralsOfForm ) {

	//  	edgeIntegrals.getCol( indexOfForm, tmpEdgeIntegralVector );

	//  	cb.computeBPeriodsForEdgeIntegrals( edgeIntegralsOfForm, tmpEdgeIntegralVector );
	//      }

	//      public static void symmetrizePeriodMatrix( ComplexMatrix periodMatrix ) {

	//  	int genus = periodMatrix.getNumRows();

	//  	for( int i=0; i<genus; i++ )
	//  	    for( int j=i+1; j<genus; j++ ) {
	//  		periodMatrix.re[i][j] = ( periodMatrix.re[i][j] + periodMatrix.re[j][i] )/2;
	//  		periodMatrix.im[j][i] = ( periodMatrix.im[i][j] + periodMatrix.im[j][i] )/2;

	//  		periodMatrix.re[j][i] =   periodMatrix.re[i][j];
	//  		periodMatrix.im[j][i] =   periodMatrix.im[i][j];

	//  	    }
	//      }

	//      public ComplexMatrix getAPeriods() {
	//  	return aPeriods;
	//      }

	//      public ComplexMatrix getBPeriods() {
	//  	return bPeriods;
	//      }

	//      public ComplexMatrix getFormTransform() {
	//  	return formTransform;
	//      }

	//      public ComplexMatrix getPeriodMatrix() {
	//  	return periodMatrix;
	//      }

	//      void updateMatrices() {

	//  	GuIntegerMatrix.matrixMult( cb.getInitialACycles(), edgeIntegrals, aPeriods, cb.getNumOfEdges(), genus, genus );
	//  	GuIntegerMatrix.matrixMult( cb.getInitialBCycles(), edgeIntegrals, bPeriods, cb.getNumOfEdges(), genus, genus );

	//  	aPeriods.assignTranspose();
	//  	bPeriods.assignTranspose();

	//  	formTransform.assignInvert( aPeriods );
	//  	formTransform.assignTimes( new Field.Complex( 0, 2*Math.PI ) );

	//  	periodMatrix.assignTimes( formTransform, bPeriods );

	//  	if( symmetrizePeriodMatrix )
	//  	    symmetrizePeriodMatrix( periodMatrix );
	//      }
}
