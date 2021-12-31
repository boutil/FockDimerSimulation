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

package de.jtem.riemann.analyticContinuation;

import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;

public class AnalyticContinuousIntegrator extends AnalyticIntegrator {


	final AnalyticContinuousIntegratorODE ode;
	
	int numberOfFunctions = 0;
	
	public AnalyticContinuousIntegrator( AlgebraicCurve curve ) {
		
		super(curve, new AnalyticContinuousIntegratorODE( curve ) );
		
		ode = (AnalyticContinuousIntegratorODE)(super.ode);
	}
	
	void setFunctionSet( FunctionSetOnAlgebraicCurve f ) {
		ode.setFunctions(f);
		
		setNumberOfFunctions( f.getNumberOfFunctions() );
	}
	
	void setNumberOfFunctions(int numberOfFunctions) {
		if( this.numberOfFunctions == numberOfFunctions ) 
			return;
		
		this.numberOfFunctions = numberOfFunctions;
		integrator.setNumOfEquations( numberOfFunctions * 2 );
		I = new double[5 + numberOfFunctions * 2];
	}

	void startContinuation( Complex lambda, Complex mu ) {
		super.startContinuation( lambda, mu );
		
		for (int i = 0, j = 5; i < numberOfFunctions; i++) {
			I[j++] = 0.0;
			I[j++] = 0.0;
		}
	}

	void startIntegration( FunctionSetOnAlgebraicCurve f, Complex lambda, Complex mu ) {

		setFunctionSet( f );
		
		startContinuation(lambda, mu );
	}

	/**
	 * 
	 * @param f set of functions to integrate
	 * @param lambdaAtStartPoint lambda coordinate of start point (input)
	 * @param muAtStartPoint mu coordinate of start point (input)
	 * @param lambdaAtEndPoint lambda coordinate of end point (input)
	 * @param muAtEndPoint mu coordinate of end point (output)
	 * @param values integrals along presribed segment (output)
	 */
	public void integrateSegmentInLambdaPlaneTo(
			FunctionSetOnAlgebraicCurve f, 
			Complex lambdaAtStartPoint,
			Complex muAtStartPoint,
			Complex lambdaAtEndPoint,
			Complex muAtEndPoint,
			ComplexVector values ) {
		
		startIntegration( f, lambdaAtStartPoint, muAtStartPoint );
		
		continueAlongSegment( lambdaAtStartPoint, lambdaAtEndPoint, true );
		
		muAtEndPoint.assign(I[3], I[4]);
		
		writeValues( values );
	}


	/**
	 * 
	 * @param f set of functions to integrate
	 * @param lambdaAtStartPoint lambda coordinate of start point (input)
	 * @param muAtStartPoint mu coordinate of start point (input)
	 * @param lambdaAtEndPoint lambda coordinate of end point (output)
	 * @param muAtEndPoint mu coordinate of end point (input)
	 * @param values integrals along presribed segment (output)
	 */
	public void integrateSegmentInMuPlaneTo(
			FunctionSetOnAlgebraicCurve f, 
			Complex lambdaAtStartPoint,
			Complex muAtStartPoint,
			Complex lambdaAtEndPoint,
			Complex muAtEndPoint,
			ComplexVector values ) {
		
		startIntegration( f, lambdaAtStartPoint, muAtStartPoint );
		
		continueAlongSegment( muAtStartPoint, muAtEndPoint, false );
		
		lambdaAtEndPoint.assign(I[1], I[2]);
		
		writeValues( values );
	}

	
	public final void integrateAlongPathOnLambdaPlane(
			FunctionSetOnAlgebraicCurve f,
			ComplexVector path,
		Complex muAtStartPoint,
		Complex muAtEndPoint, ComplexVector values ) {
		integrateAlongPath( f, path, muAtStartPoint, muAtEndPoint, true, values );
	}
	

	public final void integrateAlongPathOnMuPlane(
			FunctionSetOnAlgebraicCurve f,
			ComplexVector path,
			Complex lambdaAtStartPoint,
			Complex lambdaAtEndPoint, ComplexVector values ) {
		integrateAlongPath( f, path, lambdaAtStartPoint, lambdaAtEndPoint, false, values );
	}
	
	final void integrateAlongPath(
			FunctionSetOnAlgebraicCurve f,
				ComplexVector path,
				Complex startPoint,
				Complex endPoint, 
				final boolean continuteOnLambdaPlane, ComplexVector values ) {
		
		setFunctionSet( f );
		
		super.continueAlongPath(path, startPoint, endPoint, continuteOnLambdaPlane );
		
		writeValues( values );
	}


	private final void writeValues(ComplexVector aVec) {

			aVec.newSize(numberOfFunctions);

			for (int i = 0, j = 5; i < numberOfFunctions; i++) {
				aVec.set(i, I[j++], I[j++]);
			}
		}


}
