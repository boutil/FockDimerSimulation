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
 * Example: An Implementation of an Abelian function.
 * <p>
 * Let <code>a,b,c,d <small>&isin;</small> C<sup>g</sup></code> with 
 * <code>a+b=c+d mod 2&pi;iZ<sup>g</sup></code>, then 
 * <p align=center>
 * <code>
 * <table>
 * <tr> <td>              </td> <td> &theta;(z+a|B)&middot;&theta;(z+b|B) </td> </tr>
 * <tr> <td>f(z) = &nbsp; </td> <td> <hr noshade size=1>                  </td> </tr>
 * <tr> <td>              </td> <td> &theta;(z+c|B)&middot;&theta;(z+d|B) </td> </tr>
 * </table>
 * </code>
 * <p>
 * defines an Abelian function, which is represented by this class. 
 * To evaluate it you only need to envoke {@link #applyTo(ComplexVector)} or
 * {@link #applyTo(ComplexVector,Field.Complex)}. Below the implementation of the last:
 * <pre>
 *   public void applyTo( ComplexVector z, Field.Complex result ) {
 * 
 *	zPlusA.assignPlus( z, a );
 *	zPlusB.assignPlus( z, b );
 *	zPlusC.assignPlus( z, c );
 *	zPlusD.assignPlus( z, d );
 *
 *	theta.theta( zPlusA, factorOfZPlusA, thetaSumOfZPlusA ); 
 *	theta.theta( zPlusB, factorOfZPlusB, thetaSumOfZPlusB ); 
 *	theta.theta( zPlusC, factorOfZPlusC, thetaSumOfZPlusC ); 
 *	theta.theta( zPlusD, factorOfZPlusD, thetaSumOfZPlusD ); 
 *                                 	
 *	factor.assign(      factorOfZPlusA );                       
 *	factor.assignPlus(  factorOfZPlusB );                       
 *	factor.assignMinus( factorOfZPlusC );                       
 *	factor.assignMinus( factorOfZPlusD );                       
 *                                                            
 *	expOfFactor.assignExp( factor );
 *
 *	result.assign(       expOfFactor      );                                 
 *	result.assignTimes(  thetaSumOfZPlusA );                    
 *	result.assignTimes(  thetaSumOfZPlusB );
 *	result.assignDivide( thetaSumOfZPlusC );
 *	result.assignDivide( thetaSumOfZPlusD );
 *   }
 * </pre>
 * <p>
 * <code>AbelianFunction</code> is executable. Its {@link #main(String[])} performs
 * a periodicity test by envoking
 * {@link #testPeriodicity(ComplexVector,boolean)} on a 100 different random arguments.
 * @see Theta
 * @author Markus Schmies
 */
public class AbelianFunction {

    Theta theta;

    ComplexVector a, b, c, d;

    /**
     * Creates Abelean function with prescribed <code>periodMatrix</code>,
     * vectors <code>a,b,c,d</code> and, error tolerance <code>tol</code>.
     * The vectors <code>a,b,c,d <small>e</small> C<sup>g</sup></code> have to fill:
     * <p align=center>
     * <code>
     * a+b=c+d mod 2&pi;iZ<sup>g</sup>
     * </code>
     * Unfortunately <code>tol</code> is only the error tolerance for the four
     * individual theta sums. Use {@link #getErrorOfLastEvaluation} to recieve
     * an error estimate for the Abelian function, which depends on the argument.
     * See the description of {@link Theta} for a more detailed discussion.
     * @param periodMatrix symmetric complex matrix with negative definite real part
     * @param a vector compatible with <code>periodMatrix</code> 
     * @param b vector compatible with <code>periodMatrix</code> 
     * @param c vector compatible with <code>periodMatrix</code> 
     * @param d vector compatible with <code>periodMatrix</code> 
     * @param tol positive number
     */
    public AbelianFunction( ComplexMatrix periodMatrix,
			    ComplexVector a, 
			    ComplexVector b, 
			    ComplexVector c, 
			    ComplexVector d, double tol ) {

	theta = new Theta( periodMatrix, tol );

	this.a = new ComplexVector( a );
	this.b = new ComplexVector( b );
	this.c = new ComplexVector( c );
	this.d = new ComplexVector( d );

	ComplexVector res = a.plus(b).minus(c).minus(d).times( new Complex( 0, 2*Math.PI ) );

	//if( !res.round().equals( res ) ) 
	//  throw new IllegalArgumentException( "vectors do not fulfill criteria" );
    }


    /**
     * Applies Abelian function to <code>z</code>.
     * @param z agrument vector
     * @return value of Abelian function at <code>z</code>
     */
    public Complex applyTo( ComplexVector z ) {
	Complex result = new Complex();
	applyTo( z, result );
	return result;
    }

    Complex factor      = new Complex();
    Complex expOfFactor = new Complex();
    
    Complex factorOfZPlusA = new Complex();
    Complex factorOfZPlusB = new Complex();
    Complex factorOfZPlusC = new Complex();
    Complex factorOfZPlusD = new Complex();
    
    Complex thetaSumOfZPlusA = new Complex();
    Complex thetaSumOfZPlusB = new Complex();
    Complex thetaSumOfZPlusC = new Complex();
    Complex thetaSumOfZPlusD = new Complex();

    ComplexVector zPlusA = new ComplexVector();
    ComplexVector zPlusB = new ComplexVector();
    ComplexVector zPlusC = new ComplexVector();
    ComplexVector zPlusD = new ComplexVector();

    /**
     * Applies Abelian function to <code>z</code>.
     * @param z agrument vector
     * @param result value of Abelian function at <code>z</code>
     */
    public void applyTo( ComplexVector z, Complex result ) {

	zPlusA.assignPlus( z, a );
	zPlusB.assignPlus( z, b );
	zPlusC.assignPlus( z, c );
	zPlusD.assignPlus( z, d );

	theta.theta( zPlusA, factorOfZPlusA, thetaSumOfZPlusA ); 
	theta.theta( zPlusB, factorOfZPlusB, thetaSumOfZPlusB ); 
	theta.theta( zPlusC, factorOfZPlusC, thetaSumOfZPlusC ); 
	theta.theta( zPlusD, factorOfZPlusD, thetaSumOfZPlusD ); 
                                 	
	factor.assign(      factorOfZPlusA );                       
	factor.assignPlus(  factorOfZPlusB );                       
	factor.assignMinus( factorOfZPlusC );                       
	factor.assignMinus( factorOfZPlusD );                       
                                                            
	expOfFactor.assignExp( factor );

	result.assign(       expOfFactor      );                                 
	result.assignTimes(  thetaSumOfZPlusA );                    
	result.assignTimes(  thetaSumOfZPlusB );
	result.assignDivide( thetaSumOfZPlusC );
	result.assignDivide( thetaSumOfZPlusD );
    }

    /**
     * Returns error of last Evaluation.
     */
    public double getErrorOfLastEvaluation() {
	return expOfFactor.abs() * theta.getTol() *
	    ( thetaSumOfZPlusA.abs() + thetaSumOfZPlusA.abs() ) /
	    ( thetaSumOfZPlusC.abs() * thetaSumOfZPlusD.abs() );
    }

    /**
     * Performs periodicity test.
     * The test is quite simple. It adds each of the <code>2g</code> vectors of the basis
     * of the lattice to the argument <code>Z</code> and checks if all these
     * results coincide.
     * @param Z argument vector
     * @param verbose controles wheter the check is done verbose or silent.
     */
    public void testPeriodicity( final ComplexVector Z, final boolean verbose ) {
	
	ComplexVector z = new ComplexVector( Z );
	
	Complex valueAtZ = applyTo( Z );

	ComplexMatrix B = theta.getPeriodMatrix();

	int dim = theta.getDim();

	if( verbose ) {
	    z.print( "test periodicity at: " );
	    System.out.println( "value is:       "  + valueAtZ );
	}

	for( int i=0; i<dim; i++ ) {
	    z.setIm( i, z.getIm( i ) + 2 * Math.PI );
	    if( verbose )
		System.out.println( "value at shift: " + applyTo( z ) );
	    if( valueAtZ.dist( applyTo( z ) ) > getErrorOfLastEvaluation() )
		throw new RuntimeException( "periodicity test failed" );
	}

	for( int i=0; i<dim; i++ ) {
	    z.assignPlus( Z, B.getCol(i) );
	    if( verbose )
		System.out.println( "value at shift: " + applyTo( z ) );
	    if( valueAtZ.dist( applyTo( z ) ) > getErrorOfLastEvaluation() )
		throw new RuntimeException( "periodicity test failed" );
	}
    }

    /**
     * Performs a periodicity test for a 2-dimensional example by envoking
     * {@link #testPeriodicity(ComplexVector,boolean)} on a 100 different random arguments.
     * @param argv has no effect
     */
    public static void main( String [] argv ) { 

	ComplexMatrix B = new ComplexMatrix( 2 ); 

	B.set( 0,0, 1.690983006, 0.9510565162 ); 
	B.set( 0,1, 1.500000000, 0.3632712640 );

	B.set( 1,0, 1.500000000, 0.3632712640 ); 
	B.set( 1,1, 1.309016994, 0.9510565162 );
	
	B.assignTimes( new Complex( 0, 2 * Math.PI ) );

	ComplexVector a = new ComplexVector( 2 );
	ComplexVector b = new ComplexVector( 2 );
	ComplexVector c = new ComplexVector( 2 );
	ComplexVector d = new ComplexVector( 2 );

	a.assignRandom();
	b.assignRandom();
	c.assignRandom();
	d.assignMinus( a.plus( b ), c );

	d.im[0] += 2*Math.PI;

	AbelianFunction abelian = new AbelianFunction( B, a, b, c, d, 1e-10);

	ComplexVector z = new ComplexVector( 2 );
	
	try {
	    for( int i=0; i<100; i++ ) {
		z.assignRandom();
		abelian.testPeriodicity( z, true );
	    }
	    System.out.println( "Test was successful!" );
	}
	catch( RuntimeException e ) {
	    System.out.println( "test failed" );
	}
    }

}




