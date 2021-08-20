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

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.constrainedComplex.Fix;
import de.jtem.riemann.constrainedComplex.UnitCircleSymmetrie;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.surface.Origin;

public class HyperEllipticCurveExample extends HyperEllipticCurve {

    private static final long serialVersionUID = 1L;
    
    public HyperEllipticCurveExample() {

	branchPoint = new BranchPoint[ 5 ];
	
	Complex p = new Complex( 0.1413, 0.1018 );
	Complex q = p.invert().conjugate();

	branchPoint[0] = new BranchPoint( 0, 0 );

	new Fix( branchPoint[0] );

	branchPoint[1] = new BranchPoint( p.re, p.im );
	branchPoint[2] = new BranchPoint( q.re, q.im );
	branchPoint[3] = new BranchPoint( p.re,-p.im );
	branchPoint[4] = new BranchPoint( q.re,-q.im );

	new UnitCircleSymmetrie( branchPoint[1], branchPoint[2] );
	new UnitCircleSymmetrie( branchPoint[3], branchPoint[4] );

	origin = new Origin( 0.5, 0 );

	update();
    }	

}









    





