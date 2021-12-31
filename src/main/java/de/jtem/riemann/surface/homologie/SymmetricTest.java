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

package de.jtem.riemann.surface.homologie;

import de.jtem.blas.IntegerMatrix;
import de.jtem.blas.IntegerVector;

class SymmetricTest {

    public static void main(String[] args) throws java.io.IOException {

	if( args.length > 0 && args[0].equals("real") )
	    testRealLocus();
	else
	    testBasis();
	
    }

    static void testBasis() {
	IntegerMatrix monodr = new IntegerMatrix ( 9, 4 );
	
	monodr.set( 0, 0, 2 );
	monodr.set( 0, 1, 1 );
	monodr.set( 0, 2, 0 );
	monodr.set( 0, 3, 3 );
	
	monodr.set( 1, 0, 3 );
	monodr.set( 1, 1, 2 );
	monodr.set( 1, 2, 0 );
	monodr.set( 1, 3, 1 );
	
	monodr.set( 2, 0, 2 );
	monodr.set( 2, 1, 1 );
	monodr.set( 2, 2, 0 );
	monodr.set( 2, 3, 3 );
	
	monodr.set( 3, 0, 1 );
	monodr.set( 3, 1, 2 );
	monodr.set( 3, 2, 0 );
	monodr.set( 3, 3, 3 );

	monodr.set( 4, 0, 1 );
	monodr.set( 4, 1, 0 );
	monodr.set( 4, 2, 3 );
	monodr.set( 4, 3, 2 );

	monodr.set( 5, 0, 1 );
	monodr.set( 5, 1, 3 );
	monodr.set( 5, 2, 2 );
	monodr.set( 5, 3, 0 );

	monodr.set( 6, 0, 3 );
	monodr.set( 6, 1, 1 );
	monodr.set( 6, 2, 2 );
	monodr.set( 6, 3, 0 );

	monodr.set( 7, 0, 3 );
	monodr.set( 7, 1, 0 );
	monodr.set( 7, 2, 2 );
	monodr.set( 7, 3, 1 );

	monodr.set( 8, 0, 3 );
	monodr.set( 8, 1, 1 );
	monodr.set( 8, 2, 2 );

	CanonicalBasis transf = new CanonicalBasis( monodr, null, null );

	IntegerVector tau0 = new IntegerVector (4);
	tau0.set( 0, 0 );
	tau0.set( 1, 1 );
	tau0.set( 2, 3 );
	tau0.set( 3, 2 );
	SymmetricBasis.test( transf , 2 , 3 , 2, tau0 );
	
    }
    
    static void testRealLocus() {
	
	IntegerMatrix monodr = new IntegerMatrix ( 9, 4 );
	
	monodr.set( 0, 0, 3 );
	monodr.set( 0, 1, 2 );
	monodr.set( 0, 2, 1 );
	monodr.set( 0, 3, 0 );
	
	monodr.set( 1, 0, 3 );
	monodr.set( 1, 1, 2 );
	monodr.set( 1, 2, 1 );
	monodr.set( 1, 3, 0 );
	
	monodr.set( 2, 0, 2 );
	monodr.set( 2, 1, 1 );
	monodr.set( 2, 2, 0 );
	monodr.set( 2, 3, 3 );
	
	monodr.set( 3, 0, 1 );
	monodr.set( 3, 1, 2 );
	monodr.set( 3, 2, 0 );
	monodr.set( 3, 3, 3 );
	
	monodr.set( 4, 0, 2 );
	monodr.set( 4, 1, 3 );
	monodr.set( 4, 2, 0 );
	monodr.set( 4, 3, 1 );
	
	monodr.set( 5, 0, 1 );
	monodr.set( 5, 1, 0 );
	monodr.set( 5, 2, 3 );
	monodr.set( 5, 3, 2 );
	
	monodr.set( 6, 0, 1 );
	monodr.set( 6, 1, 3 );
	monodr.set( 6, 2, 0 );
	monodr.set( 6, 3, 2 );
	
	monodr.set( 7, 0, 1 );
	monodr.set( 7, 1, 3 );
	monodr.set( 7, 2, 2 );
	monodr.set( 7, 3, 0 );
	
	monodr.set( 8, 0, 0 );
	monodr.set( 8, 1, 3 );
	monodr.set( 8, 2, 2 );
	monodr.set( 8, 3, 1 );
	
	CanonicalBasis transf = new CanonicalBasis( monodr, null, null );
	    
	IntegerVector tau0 = new IntegerVector (4);
	tau0.set( 0, 1 );
	tau0.set( 1, 0 );
	tau0.set( 2, 3 );
	tau0.set( 3, 2 );
	SymmetricBasis.test( transf , 2 , 3 , 2, tau0 );
	
    }
}







