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

import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.calculus.odeSolving.ODE;
import de.jtem.numericalMethods.calculus.odeSolving.ODEIntermediateResultListener;

public class AnalyticContinuationODE implements ODE {
	
	int numOfEvaluations = 0;

	final Corrector cf;
	final AlgebraicCurve curve;
	
	boolean continueInLambdaPlane;
	boolean continueInMuPlane;   //TODO: remove redundant field

	public AnalyticContinuationODE( AlgebraicCurve curve ) {
		this( curve, true );
	}
	
	public AnalyticContinuationODE( AlgebraicCurve curve, boolean continueInLambdaPlane ) {
		 this.curve = curve;
		 	 
		 cf = new Corrector(curve);
		 
		 setContinueInLambdaPlane( continueInLambdaPlane );
	}
	
	void setContinueInLambdaPlane( boolean continueInLambdaPlane ) {
		this.continueInLambdaPlane = continueInLambdaPlane;
		this.continueInMuPlane = !continueInLambdaPlane;
	}
	
	final ODEIntermediateResultListener fcf = new ODEIntermediateResultListener() {

		private static final long serialVersionUID = 1L;

		final Complex lambda = new Complex();
		final Complex mu = new Complex();


		public void intermediateResult(double t, double[] x, int offsetX) {

			lambda.assign(x[0 + offsetX], x[1 + offsetX]);
			mu.assign(x[2 + offsetX], x[3 + offsetX]);

			if (continueInMuPlane)
				cf.correctLambda(lambda, mu);
			else
				cf.correctMu(lambda, mu);

			x[0 + offsetX] = lambda.re;
			x[1 + offsetX] = lambda.im;

			x[2 + offsetX] = mu.re;
			x[3 + offsetX] = mu.im;
		}

	};

	public int getNumberOfEquations() {
		return 4;
	}


	final Complex lambda = new Complex();
	final Complex mu = new Complex();

	final Complex direction = new Complex();
	final Complex last = new Complex();

	final Complex Q = new Complex();
	final Complex dQdLambda = new Complex();
	final Complex dQdMu = new Complex();

	final Complex dLambdaDt = new Complex();
	final Complex dMuDt = new Complex();


	double localRadiusSqrOfLast = 0;

	public void eval(
		double t,
		double[] x,
		double[] y) {

		lambda.assign(x[0], x[1]);
		mu.assign(x[2], x[3]);

		if (continueInMuPlane) {

			cf.correctLambda(lambda, mu);

			if (false && last.distSqr(lambda) > localRadiusSqrOfLast) {
				System.out.println(
					"localRadiusSqrOfQ = " + localRadiusSqrOfLast);
				System.out.println("dist = " + last.distSqr(lambda));
				System.out.println("lambda = " + lambda);
				System.out.println("branch point = " + last);
				System.out.println("t = " + t);

				// this is a big problem
				throw new RuntimeException("exit local region");
			}

			curve.evaluateCurve(lambda, mu, Q, dQdLambda, dQdMu);

			dMuDt.assign(direction);

			dLambdaDt.assignDivide(dQdMu, dQdLambda);
			dLambdaDt.assignTimes(-1);
			dLambdaDt.assignTimes(dMuDt);

		} else {

			cf.correctMu(lambda, mu);

			curve.evaluateCurve(lambda, mu, Q, dQdLambda, dQdMu);

			dLambdaDt.assign(direction);

			dMuDt.assignDivide(dQdLambda, dQdMu);
			dMuDt.assignTimes(-1);
			dMuDt.assignTimes(dLambdaDt);

		}

		if (continueInMuPlane) {
			y[0] = dLambdaDt.re;
			y[1] = dLambdaDt.im;

			y[2] = direction.re;
			y[3] = direction.im;
		} else {
			y[0] = direction.re;
			y[1] = direction.im;

			y[2] = dMuDt.re;
			y[3] = dMuDt.im;
		}

		numOfEvaluations++;
		
	}

	
}
