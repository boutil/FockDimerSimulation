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

import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.calculus.function.RealFunctionOfSeveralVariables;
import de.jtem.numericalMethods.calculus.minimizing.Info;
import de.jtem.numericalMethods.calculus.minimizing.NelderMead;
import de.jtem.riemann.constrainedComplex.IllegalCoordinateException;

public abstract class RamifiedCoveringVariation extends RamifiedCovering
    implements RealFunctionOfSeveralVariables {

    private static final long serialVersionUID = 1L;

    boolean returnBigDoubleValue;

    boolean weAreMinimizing = false;

    protected double  doubleValue = 123456789, minDoubleValue;

    RamifiedCoveringState minDoubleValueState;

    private final Complex tmp = new Complex();

    public RamifiedCoveringState getMinDoubleValueState() {
	return minDoubleValueState;
    }

    public double getDoubleValue() {

	if( returnBigDoubleValue )
	    return 1e100;

	try {
	    update();
	} catch( RuntimeException e ) {

	    System.out.println( "caught exception: " + e.getMessage() );

	    double [] value = getDoubleArrayValue();

	    for( int i=0; i<value.length; i++ )
		System.out.print( value[i] + "  " );

	    System.out.println( "end of uniformization data");

  	    setState( minDoubleValueState );

	    //throw e;

	    returnBigDoubleValue = true;
	    return 1e100;
	}

	if( weAreMinimizing && doubleValue < minDoubleValue ) {

	    getState( minDoubleValueState );

	    //firePropertyChange( null, true );

	    minDoubleValue = doubleValue;

	    test();

	    //System.out.println( "minDoubleValue = " + minDoubleValue );
	}

	return doubleValue;
    }

    public int getDoubleArrayParameterLength() {
      return getNumberOfVariables();
    }

    public void setByParameter( double [] p ) {
	setDoubleArrayParameter ( p, 0 );
    }

    public int getNumberOfVariables() {
      int numOfParameters = 0;
        for( int i=0; i<numOfSurfacePoints; i++ )
            numOfParameters += surfacePoint[i].getNumOfParameters();

        return numOfParameters;
    }

    public double eval( double [] data ) {
      setDoubleArrayParameter( data, 0 );
      return getDoubleValue();
    }

    double factor = 100;

    public double getFactor() {
	return factor;
    }

    public void setFactor( double aFactor ) {
	factor = aFactor;
    }

    public void setDoubleArrayParameter( double [] p, int offset ) {

	returnBigDoubleValue = false;

	RamifiedCoveringState saveState;

	if( weAreMinimizing ) {
	    saveState = minDoubleValueState;

	    setState( minDoubleValueState );

	} else
	    saveState = getState();

	try {

	    for( int i=0; i<numOfBranchPoints; i++ ) {
		int newOffset = offset +  branchPoint[i].getNumOfParameters();
		if( offset < newOffset )
		    branchPoint[i].setByParameter( p, offset );
		offset = newOffset;
	    }

	    for( int i=0; i<numOfSingularPoints; i++ ) {
		int newOffset = offset +  singularPoint[i].getNumOfParameters();
		if( offset < newOffset )
		    singularPoint[i].setByParameter( p, offset );
		offset = newOffset;
	    }

	    for( int i=0; i<numOfDistinguishedPoints; i++ ) {
		int newOffset = offset +  distinguishedPoint[i].getNumOfParameters();
		if( offset < newOffset )
		    distinguishedPoint[i].setByParameter( p, offset );
		offset = newOffset;
	    }

	    outdate();

	} catch( IllegalCoordinateException e ) {

	    setState( saveState );

	    if( weAreMinimizing ) {
		returnBigDoubleValue = true;

		outdate();
	    }
	}

    }


    public double []  getDoubleArrayValue() {
	double [] value = new double[ getDoubleArrayParameterLength() ];

	getValue( value );

	return value;
    }

    public int getDoubleArrayValueLength() {
	return getDoubleArrayParameterLength();
    }

    public void getValue( double [] value ) {
	getDoubleArrayValue ( value, 0 );
    }

    public void getDoubleArrayValue( double [] value, int offset ) {


	    for( int i=0; i<numOfBranchPoints; i++ ) {
		int newOffset = offset +  branchPoint[i].getNumOfParameters();
		if( offset < newOffset )
		    branchPoint[i].getValue( value, offset );
		offset = newOffset;
	    }

	    for( int i=0; i<numOfSingularPoints; i++ ) {
		int newOffset = offset +  singularPoint[i].getNumOfParameters();
		if( offset < newOffset )
		    singularPoint[i].getValue( value, offset );
		offset = newOffset;
	    }

	    for( int i=0; i<numOfDistinguishedPoints; i++ ) {
		int newOffset = offset +  distinguishedPoint[i].getNumOfParameters();
		if( offset < newOffset )
		    distinguishedPoint[i].getValue( value, offset );
		offset = newOffset;
	    }

    }

    public void test() {

    }

    public void minimize() {
	minimize( 50, 1e-6 );
    }

    public void minimize( int n, double ftol ) {

	weAreMinimizing = true;

	try {


	    setEnablePropertyChange( false );

	    update();

	    double [] value = getDoubleArrayValue();

	    minDoubleValue      = getDoubleValue();
	    minDoubleValueState = getState();

	    //Powell.search( value, ftol, n, this, new Info(true) );
	    NelderMead.search( value, ftol, n, this, new Info(true) );
	    setState( minDoubleValueState );

	    setEnablePropertyChange( true );

	    update();

	    System.out.println( "minDoubleValue = " + minDoubleValue );

	    test();

	} finally {
	    returnBigDoubleValue = weAreMinimizing = false;
	}
    }
}














