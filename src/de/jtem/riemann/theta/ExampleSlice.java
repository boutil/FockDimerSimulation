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

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;

/**
 * Example from the paper <em>Computing Riemann Theta Functions</em>.
 * <p>
 * The example computes the following slice on a <code>101x101</code> sample:
 * <p align=center>
 *   <code>
 *     (x,y) |-> thetaSum( (x+iy,0) |B)
 *   </code>.
 * </p>
 * Attention! This example has no visual output.
 * </p>
 * Below the source code of this example:
 * </p>
 * <pre>
 *  package de.jtem.riemann.theta;
 * 
 *  import de.jtem.mfc.field.Field;
 *  import de.jtem.blas.ComplexVector;
 *  import de.jtem.blas.ComplexMatrix;
 *  
 *  public class ExampleSlice {
 * 
 *     public static void main( String [] argv ) { 
 * 
 *        Field.Complex factor = new Field.Complex();
 *        Field.Complex sum    = new Field.Complex();
 *
 *        // create the period matrix
 *        ComplexMatrix B = new ComplexMatrix( 2 ); 
 * 
 *        // set the period matrix
 * 
 *        B.set( 0,0, 1.690983006, 0.9510565162 ); 
 *        B.set( 0,1, 1.500000000, 0.3632712640 );
 * 
 *        B.set( 1,0, 1.500000000, 0.3632712640 ); 
 *        B.set( 1,1, 1.309016994, 0.9510565162 );
 * 
 *        // use the different normalization
 *        B.assignTimes( new Field.Complex( 0, 2 * Math.PI ) );
 * 
 *        // create a complex argument vector
 *        ComplexVector V = new ComplexVector( 2 );
 * 
 *        // create a theta function instance with period matrix 
 *        Theta theta = new Theta( B, 0.0001 ); 
 * 
 *        // create storage for the slice
 *        Field.Complex[][] slice = new Field.Complex[101][101]; 
 * 
 *        // create dummy variable for the factor
 *        Field.Complex factor = new Field.Complex();
 *
 *        // loop over grid
 *        for( int i=0; i<101; i++ )
 * 
 *        for( int j=0; j<101; j++ ) {
 * 
 * 	     // set grid vector
 * 	     V.set( 0, 2 * i * Math.PI / 100, 2 * j *  Math.PI / 100 ); 
 * 
 * 	     // evaluate theta function at grid vector
 * 	     theta.theta( V, factor, slice[i][j] ); 
 *        }
 *     }
 *  }
 * </pre>
 * @see Theta
 * @author Markus Schmies
 */
public class ExampleSlice {

    public static void main( String [] argv ) { 

	// create the period matrix
	ComplexMatrix B = new ComplexMatrix( 2 ); 

	// set the period matrix
	B.set( 0,0, 1.690983006, 0.9510565162 );
	B.set( 0,1, 1.500000000, 0.3632712640 );

	B.set( 1,0, 1.500000000, 0.3632712640 ); 
	B.set( 1,1, 1.309016994, 0.9510565162 );

	// use the different normalization
	B.assignTimes( new Complex( 0, 2 * Math.PI ) );

	// create a complex argument vector
	ComplexVector V = new ComplexVector( 2 );

	// create a theta function instance with period matrix 
	Theta theta = new Theta( B, 0.0001 ); 

	// create storage for the slice
	Complex[][] slice = new Complex[101][101]; 

	// create variables for evaluating the theta functin
	Complex factor   = new Complex();	
	
	// loop over grid
	for( int i=0; i<101; i++ )
	    for( int j=0; j<101; j++ ) {

		// set grid vector
		V.set( 1, -10 * i * Math.PI / 100, 2 * j * Math.PI / 100 ); 

		// evaluate theta function at grid vector
		theta.theta( V, factor, slice[i][j] ); 
	    }
    }

}




