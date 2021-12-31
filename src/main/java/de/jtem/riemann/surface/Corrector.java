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

import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson;
import de.jtem.riemann.analyticContinuation.AlgebraicCurve;

class Corrector implements NewtonRaphson.ComplexFunctionWithDerivative, Serializable, Cloneable {
	

	/**
	 * 
	 */
	private final AlgebraicCurve curve;

	/**
	 * @param curve
	 */
	Corrector(AlgebraicCurve curve) {
		this.curve = curve;
	}



	private static final long serialVersionUID = 1L;

	boolean correctMu = true;

	final Complex lambda = new Complex();
	final Complex mu = new Complex();

	final Complex Q = new Complex();
	final Complex dQdLambda = new Complex();
	final Complex dQdMu = new Complex();

	final double[] rootValues = new double[6];

	void correctMu(final Complex aLambda, final Complex aMu) {
		correctMu = true;
		lambda.assign(aLambda);

		try {
			de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.search( this, aMu.re, aMu.im, rootValues);
		} catch (RuntimeException e) {
			de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.DEBUG = true;

			try {
				de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.search( this, aMu.re, aMu.im, rootValues);
			} catch (RuntimeException ee) {
			}

			de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.DEBUG = false;
			System.out.println("could not solve equation");

			throw e;
		}

		aMu.assign(rootValues[0], rootValues[1]);
	}

	void correctLambda(final Complex aLambda, final Complex aMu) {
		correctMu = false;
		mu.assign(aMu);

		try {
			de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.search(this, aLambda.re, aLambda.im, rootValues);
		} catch (RuntimeException e) {
			de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.DEBUG = true;

			try {
				de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.search(this, aLambda.re, aLambda.im, rootValues);
			} catch (RuntimeException ee) {
			}

			de.jtem.numericalMethods.calculus.rootFinding.NewtonRaphson.DEBUG = false;
			System.out.println("could not solve equation");
			System.out.println("lambda=" + lambda);
			System.out.println("mu    =" + mu);

			throw e;
		}

		aLambda.assign(rootValues[0], rootValues[1]);
	}



	public void eval(double x, double y, double[] f, int offsetF, double[] df, int offsetDF) {

		if (correctMu)
			mu.assign(x,y);
		else
			lambda.assign(x,y);

		this.curve.evaluateCurve(lambda, mu, Q, dQdLambda, dQdMu);
		

		f[offsetF+0] = Q.re;
		f[offsetF+1] = Q.im;

		if (correctMu) {
			df[offsetDF+0] = dQdMu.re;
			df[offsetDF+1] = dQdMu.im;
		} else {
			df[offsetDF+0] = dQdLambda.re;
			df[offsetDF+1] = dQdLambda.im;
		}
	}
	
}