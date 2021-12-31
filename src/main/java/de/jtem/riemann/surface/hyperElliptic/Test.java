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

package de.jtem.riemann.surface.hyperElliptic;

import de.jtem.blas.ComplexMatrix;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;

class Test {

    public static void main(String[] args) throws java.io.IOException {

	/* genus 1 example
	BranchPoint [] branchPoint = new BranchPoint[4];

	branchPoint[0] = new BranchPoint( 1,  1 );
	branchPoint[1] = new BranchPoint(-1,  1 );
	branchPoint[2] = new BranchPoint(-1, -1 );
	branchPoint[3] = new BranchPoint( 1, -1 );
	*/

	/* genus 2 example */
	BranchPoint [] branchPoint = new BranchPoint[6]; 

	branchPoint[0] = new BranchPoint( 2,  1.0 );
	branchPoint[1] = new BranchPoint( 1,  1.1 );
	branchPoint[2] = new BranchPoint( 0,  1.2 );
	branchPoint[3] = new BranchPoint( 0, -1.3 );
	branchPoint[4] = new BranchPoint( 1, -1.4 );
	branchPoint[5] = new BranchPoint( 2, -1.5 );
	
	HyperEllipticSurface aHyperEllipticSurface 
	    = new HyperEllipticSurface( branchPoint, null, null, new Complex( 0, 0 ) );

	aHyperEllipticSurface.setNumOfSteps( 301 );

	
	ComplexMatrix periodMatrix = new ComplexMatrix( aHyperEllipticSurface.getGenus() );

	aHyperEllipticSurface.computePeriodMatrix( periodMatrix );

	periodMatrix.print("period matrix");	
	
    }
    

}




