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

import de.jtem.moebiusViewer.MoebiusGraphics;
import de.jtem.moebiusViewer.shape.IndexedLineSet;


public class AbstractAlgebraicCurveShape extends RamifiedCoveringShape {

    private static final long serialVersionUID = 1L;
    
    protected IndexedLineSet indexedLineSet;
    
    /**
       * Get the value of indexedLineSet.
       * @return Value of indexedLineSet.
       */
    public IndexedLineSet getIndexedLineSet() {return indexedLineSet;}
        
    protected AbstractAlgebraicCurve simple;

    int getTotalNumberOfPoints() {

	int number = 0;

	for( int i=0;i<simple.getNumOfBranchPoints(); i++ )
	    number += simple.branchPointLoop[i].size();

	for( int i=0;i<simple.getNumOfSingularPoints(); i++ )
	    number += simple.singularPointLoop[i].size();

	for( int i=0;i<simple.getNumOfDistinguishedPoints(); i++ )
	    number += simple.distinguishedPointPath[i].size();
	
	return number;
    }

    public void computeIndexedLineSet() {

	simple = (AbstractAlgebraicCurve)surface;

	simple.updateLoops();

	int index[]   = new int[ getTotalNumberOfPoints() ];
	int length [] = new int[ surface.getNumOfSurfacePoints() ];
	double xy[]   = new double[ 2 * index.length ];

	int l = 0, k = 0;

	for( int i=0; i<index.length; i++ )
	    index[i] = i;
	
	for( int i=0;i<simple.getNumOfBranchPoints(); i++ ) {

	    length[k++] = simple.branchPointLoop[i].size();
	    
	    for( int j=0; j<simple.branchPointLoop[i].size(); j++ ) {

		xy[l++] = simple.branchPointLoop[i].re[j];
		xy[l++] = simple.branchPointLoop[i].im[j];
	    }
	}

	for( int i=0;i<simple.getNumOfSingularPoints(); i++ ) {

	    length[k++] = simple.singularPointLoop[i].size();
	    
	    for( int j=0; j<simple.singularPointLoop[i].size(); j++ ) {
		xy[l++] = simple.singularPointLoop[i].re[j];
		xy[l++] = simple.singularPointLoop[i].im[j];
	    }
	}

	for( int i=0;i<simple.getNumOfDistinguishedPoints(); i++ ) {

	    length[k++] = simple.distinguishedPointPath[i].size();
	    
	    for( int j=0; j<simple.distinguishedPointPath[i].size(); j++ ) {
		xy[l++] = simple.distinguishedPointPath[i].re[j];
		xy[l++] = simple.distinguishedPointPath[i].im[j];
	    }
	}
	
	indexedLineSet = new IndexedLineSet ( xy, index, length );
    }

    public void compute() {

	super.compute();

	computeIndexedLineSet();
    }

    public AbstractAlgebraicCurveShape() {
    }

    public AbstractAlgebraicCurveShape( RamifiedCovering aSurface ) {
	setSurface( aSurface );	
    }

    public void draw (MoebiusGraphics context) {

	if( indexedLineSet != null )
	    indexedLineSet.draw( context );

	super.draw( context );


    }  
}
 











