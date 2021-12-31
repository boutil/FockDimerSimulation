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

package de.jtem.riemann.schottky;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexConstant;


public class FixPointsOfInvolution {

    private static final Complex onePlusMu   = new Complex();
    private static final Complex oneMinusMu  = new Complex();
    private static final Complex ratio       = new Complex();
    private static final Complex root        = new Complex();
    private static final Complex ratioTimesA = new Complex();

    public static void compute( final Schottky schottky, final int anIndex,
				final Complex firstFixPoint, final Complex secondFixPoint ) {

	compute( schottky.fixpoint[anIndex][0],
		 schottky.fixpoint[anIndex][1],
		 schottky.mu[anIndex], firstFixPoint, secondFixPoint );
    }

    public static void compute( final Complex A, final Complex B, final Complex mu,
			   final Complex firstFixPoint, final Complex secondFixPoint ) {

	onePlusMu.assignPlus (  ComplexConstant.ONE, mu );
	oneMinusMu.assignMinus( ComplexConstant.ONE, mu );

	ratio.assignDivide( onePlusMu, oneMinusMu );

	root.assignTimes( ratio, ratio );
	root.assignMinus( 1 );
	root.assignTimes( A );
	root.assignTimes( A );
	root.assignSqrt();

	ratioTimesA.assignTimes( ratio, A );

	if( firstFixPoint != null )
	    firstFixPoint. assignPlus(  ratioTimesA, root );
	if( secondFixPoint != null )
	    secondFixPoint.assignMinus( ratioTimesA, root );
    }

  public static Complex [] getBranchPoints( Schottky schottky ) {
    Complex [] branchPoint = new Complex[2 * schottky.getNumGenerators() + 1];

    computeBranchPoints( schottky, branchPoint );

    return branchPoint;
  }

  public static void computeBranchPoints( Schottky schottky, Complex [] branchPoint ) {

    if( branchPoint.length != 2 * schottky.getNumGenerators() + 1 )
      throw new IllegalArgumentException( "array for branch ponits has wrong size" );

    branchPoint[0] = new Complex( 0, 0 ); // zero and infinity are always branch points

    int  numOfBranchPoints = 1;

    Complex firstFixPoint  = new Complex();
    Complex secondFixPoint = new Complex();

    Complex sqrSigmaAtZero     = new Complex();
    Complex sqrSigmaAtFixPoint = new Complex();

    Complex o = new Complex();

    for( int i=0; i<schottky.getNumGenerators(); i++ ) {

      compute( schottky, i, firstFixPoint, secondFixPoint );

      schottky.sigma( sqrSigmaAtFixPoint, firstFixPoint, o, 2 );

      System.out.println( branchPoint[ numOfBranchPoints++ ] = sqrSigmaAtFixPoint.minus( sqrSigmaAtZero ) );


      schottky.sigma( sqrSigmaAtFixPoint, secondFixPoint, o, 2 );

      System.out.println( branchPoint[ numOfBranchPoints++ ] = sqrSigmaAtFixPoint.minus( sqrSigmaAtZero ) );
    }
  }
  

  public static Complex [] getFixPoints( Schottky schottky ) {
  	Complex [] fixPoint = new Complex[2 * schottky.getNumGenerators()];

  	computeFixPoints( schottky, fixPoint );

  	return fixPoint;
  }

  public static void computeFixPoints( Schottky schottky, Complex [] fixPoint ) {

  	if( fixPoint.length != 2 * schottky.getNumGenerators()  )
  		throw new IllegalArgumentException( "array for branch ponits has wrong size" );

  	for( int i=0, j=0; i<schottky.getNumGenerators(); i++, j+=2 ) {

  		if( fixPoint[j] == null ) fixPoint[j] = new Complex();
  		if( fixPoint[j+1] == null ) fixPoint[j+1] = new Complex();
  		
  		compute( schottky, i,  fixPoint[j],  fixPoint[j+1]  );
  	}
  }
}












