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
import de.jtem.numericalMethods.calculus.odeSolving.ODE;

class AnalyticContinuousIntegratorODE extends AnalyticContinuationODE implements ODE {
		
	FunctionSetOnAlgebraicCurve functions;
	
	int numberOfFunctions;
	
	ComplexVector values = new ComplexVector();
	
	public AnalyticContinuousIntegratorODE(  AlgebraicCurve curve ) {
		super( curve );	 
	}
	
	void setFunctions( FunctionSetOnAlgebraicCurve functions ) {
		this.functions = functions;	
		
		numberOfFunctions = functions.getNumberOfFunctions();
	
		values.setSize( numberOfFunctions );
	}
	
	public int getNumberOfEquations() {
		return 4 + numberOfFunctions * 2;
	}

	public void eval(
		double t,
		double[] x,
		double[] y) {

		super.eval( t, x, y );
	
		functions.evaluateFunctions(lambda, mu, Q, dQdLambda, dQdMu, values);
	
		values.assignTimes( direction );
		for (int i = 0, j = 4; i < numberOfFunctions; i++, j+=2 ) {
			y[j   ] = values.re[i];
			y[j+1 ] = values.im[i];
		}
		
	}
}
