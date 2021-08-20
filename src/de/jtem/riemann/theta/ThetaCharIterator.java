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

import de.jtem.blas.IntegerVector;
import de.jtem.mfc.field.Complex;

/**
 * Iterates theta characteristics.
 *
 * @see ThetaWithChar
 * @author Markus Schmies
 */
public class ThetaCharIterator implements Serializable, Cloneable{

    private static final long serialVersionUID = 1L;

    ThetaWithChar theta;

    int dim = -1;

    int numOfCharacteristics;

    int alphaIndex = 0;
    int betaIndex  = 0;

    int numOfIterations;

    public ThetaCharIterator( ThetaWithChar aTheta ) {

	theta = aTheta;

	update();
    }

    /** @return 4^genus
     */
    public int getNumOfCharacteristics() {
	update();	
	return numOfCharacteristics;
    }

    /** @return 2^(genus-1) * ( 2^genus + 1 ) 
     */
    public int getNumOfEvenCharacteristics() {
	update();
	return theta.getNumOfEvenCharacteristics();
    }

    /** @return 2^(genus-1) * ( 2^genus + 1 ) 
     */
    public int getNumOfOddCharacteristics() {
	update();
	return theta.getNumOfOddCharacteristics();
    }

    
    private boolean [][] isEvenChar;

    private IntegerVector[] charVector;
    private IntegerVector y;
    private int numOfCharVectors;
    private int charVectorLength;

    private Complex tmp = new Complex();

    private void updateCharVectors() {

	charVectorLength = (int)(Math.pow( 2, dim ) + 0.5 );

	charVector = new IntegerVector[ charVectorLength ];

	y = new IntegerVector ( dim );

	numOfCharVectors = 0;

	updateCharVectors( 0 );

	isEvenChar = new boolean[numOfCharVectors][numOfCharVectors];

	for( int i=0; i<numOfCharVectors; i++ )
	    for( int j=0; j<numOfCharVectors; j++ )
		isEvenChar[i][j] = ( charVector[i].dot( charVector[j] ) ) % 2 == 0;	    
    }
    
    private void updateCharVectors( int k ) {

	if( k != dim )  {

	    y.re[k] = 0; updateCharVectors( k+1 );
	    y.re[k] = 1; updateCharVectors( k+1 );

	} else charVector[ numOfCharVectors++ ] = new IntegerVector ( y );	    
    }


    void setDim(int v) {

	if( v != dim ) {

	    dim = v;

	    updateCharVectors();

	    numOfCharacteristics = theta.getNumOfCharacteristics();

	    numOfIterations = 0;

	    alphaIndex = 0;
	    betaIndex  = 0;
	}
    }
    

    final void update() {
	theta.update();
	setDim( theta.getDim() );	
    }

    public void startCharIteration() {
	update();
	numOfIterations = 0;

	theta.alpha = charVector[ alphaIndex ];
	theta.beta  = charVector[  betaIndex ];

	theta.uptodate = false;
    }
    
    private void incrementIndices() {
	
	if( dim != theta.dim )
	    throw new RuntimeException( "during iteration genus of theta function changed" );

	betaIndex++;

	if( betaIndex == numOfCharVectors )
	    alphaIndex++;

	betaIndex  %= numOfCharVectors;
	alphaIndex %= numOfCharVectors;

	numOfIterations++;
    }

    public boolean iterateChar() {
	update();

	if( numOfIterations == numOfCharacteristics )
	    return false;

	incrementIndices();

	theta.alpha = charVector[ alphaIndex ];
	theta.beta  = charVector[  betaIndex ];

	theta.uptodate = false;

	return true;
    }

    public void startOddCharIteration() {
	update();
	numOfIterations = 0;

	while( isEvenChar[ alphaIndex ][ betaIndex ] )
	    incrementIndices();
	
	theta.alpha = charVector[ alphaIndex ];
	theta.beta  = charVector[  betaIndex ];

	theta.uptodate = false;
    }

    public boolean iterateOddChar() {
	update();
	if( numOfIterations == numOfCharacteristics )
	    return false;

	incrementIndices();

	while( numOfIterations != numOfCharacteristics && isEvenChar[ alphaIndex ][ betaIndex ] ) {
	    incrementIndices();
	}

	theta.alpha = charVector[ alphaIndex ];
	theta.beta  = charVector[  betaIndex ];

	theta.uptodate = false;

	return numOfIterations != numOfCharacteristics;
    }

    public void startEvenCharIteration() {
	update();
	numOfIterations = 0;

	while( !isEvenChar[ alphaIndex ][ betaIndex ] )
	    incrementIndices();
	
	theta.alpha = charVector[ alphaIndex ];
	theta.beta  = charVector[  betaIndex ];

	theta.uptodate = false;
    }

    public boolean iterateEvenChar() {
	update();
	if( numOfIterations == numOfCharacteristics )
	    return false;

	incrementIndices();

	while( numOfIterations != numOfCharacteristics && !isEvenChar[ alphaIndex ][ betaIndex ] ) {
	    incrementIndices();
	}

	theta.alpha = charVector[ alphaIndex ];
	theta.beta  = charVector[  betaIndex ];

	theta.uptodate = false;

	return numOfIterations != numOfCharacteristics;
    }

    private int getIndexOfCharVector( IntegerVector aCharVector ) {

	for( int i=0; i<numOfCharVectors; i++ )
	    if( aCharVector.equals( charVector[i] ) )
		return i;

	throw new IllegalArgumentException( "vector is no characteristics" );
    }

}
















