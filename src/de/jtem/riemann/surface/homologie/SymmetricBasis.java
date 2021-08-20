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

/** 
    All members of this class are defined during instantiation and stay unchanged afterwards:</br>    
    <i>H</i></br>
    <i>tau0</i></br>
    <i>realLocus</i></br>
    <i>nl</i></br>
    <i>nr</i></br>
    <i>nb</i></br>
    <br>
    During instantiation the following variables from the class <i>CanonicalBasis</i> are recalculated:</br>
    </br>
    <i>initialACycles</i></br>
    <i>initialBCycles</i></br>
    <br>
*/

public class SymmetricBasis extends CanonicalBasis implements Cloneable {

    private static final long serialVersionUID = 1L;
	
    IntegerMatrix H;

    IntegerVector tau0;

    int nl, nr, nb;

    int realLocus;

    private static boolean debug = false;

    
    public IntegerMatrix getH() {
	return new IntegerMatrix ( H );
    }

    public IntegerVector getTau0() {
	return new IntegerVector ( tau0 );
    }

    public int getNl() {
	return nl;
    }

    public int getNr() {
	return nr;
    }

    public int getNb() {
	return nb;
    }
    
    public int getRealLocus() {
	return realLocus;
    }

    public int getInvariantA() {
	return getInvariantA( H );
    }
       
    static IntegerMatrix getTransformToSymmetricBasis( IntegerMatrix tau ) {

	int numOfCycles = tau.getNumRows(), genus = numOfCycles / 2;

	IntegerMatrix F = IntegerMatrix.id( numOfCycles );
	IntegerMatrix A = F.plus( tau );
	IntegerMatrix G = new IntegerMatrix ( A );

	int zero_space = 0;

	for(int i=0; i+zero_space < numOfCycles ; i++) {

	  while( true ) {	    
	    int [] pos = new int[2];
	    int [] val = new int[2];
	    
	    GuIntegerMatrix.find2LargestInColStartingFrom( A, i, i, pos, val );
	    
	    if( pos[0] == -1 ) { 

	      zero_space++;
	      
	      if(i + zero_space >=  numOfCycles) 
		break;

	      GuIntegerMatrix.moveColToEnd ( A, i );
	      GuIntegerMatrix.moveColToEnd ( G, i );	      

	    } else {

	      if( pos[1] == -1 ) {

		if( pos[0] != i) {

		  GuIntegerMatrix.swapRows ( A, i, pos[0]);	
		  GuIntegerMatrix.swapRows ( F, i, pos[0]);	
		  
		}

		break;

	      } else {

		int k =- val[0] / val[1];

		GuIntegerMatrix.addRow ( A, pos[0] , pos[1] , k );
		GuIntegerMatrix.addRow ( F, pos[0] , pos[1] , k );
		
	      }
	    }
	  }
	}
  
	for( int i=0; i+zero_space<numOfCycles; i++) {

	  for( int j=0; j<i; j++) {

	    int k = -A.get(j,i)/A.get(i,i);

	    GuIntegerMatrix.addRow ( A, j , i , k );
	    GuIntegerMatrix.addRow ( F, j , i , k );
	    
	  }
	}

	// A.print( " A ");

	for( int i=0; i+zero_space<numOfCycles; i++) {

	  int k = A.get(i,i);

	  GuIntegerMatrix.divideCol ( A, i, k );
	  GuIntegerMatrix.divideCol ( G, i, k );     

	  for( int j=i+1; j+zero_space < numOfCycles; j++) {

	    int k1 = -A.get( i, j );

	    GuIntegerMatrix.addCol ( A, j, i, k1 );
	    GuIntegerMatrix.addCol ( G, j, i, k1 );
      
	  }              
	}
	
  
	for( int i = 0 ; i<genus; i++) {
	    for( int j=0;  j<genus; j++) {
		G.set( i,       j+genus, -F.get( j, i + genus ) );
		G.set( i+genus, j+genus,  F.get( j, i         ) );
	    }
}	

	/* End of the Step 1    */
	/* Start of the Step 2 */

	IntegerMatrix B = tau.times( G );
	B.assignPlus( G );
	A.assignTimes( F, B );

	// A.print( "  C ");  
  
	for( int i = 0 ; i<genus; i++) {
	    for( int j=0;  j<genus; j++) {
		int k = A.get( j, i + genus );
		if( k < 0  && k - 2 * (k / 2) != 0 ) k-=2;
		k=-k/2;
		GuIntegerMatrix.addCol ( G, i+genus , j , k );
	    }
	}

	/* End of the Step 2 */

	// G.print(" Tranformation matrtix ");  

	return G;
    }

    static IntegerMatrix getAntiholomorphicInvolution( IntegerMatrix tau, IntegerMatrix tr ) {

	int numOfCycles = tau.getNumRows(), genus = numOfCycles / 2;

	IntegerMatrix tr1 = new IntegerMatrix ( numOfCycles );

	for( int i = 0 ; i<genus; i++) {
	    for( int j=0;  j<genus; j++) {
		tr1.set( j,       i,        tr.get( i+genus, j+genus ) );
		tr1.set( j,       i+genus, -tr.get( i,       j+genus ) );
		tr1.set( j+genus, i,       -tr.get( i+genus, j       ) );
		tr1.set( j+genus, i+genus,  tr.get( i,       j       ) );
	    }
	}
  
	// tr1.times( tr ).print( " G^-1 * G ");
	tr1.times( tau ).assignTimes( tr );
	return tr1;
    }

    public static int getInvariantA( IntegerMatrix H ) {

	for( int i=0; i<H.getNumRows(); i++)
	    if( H.get( i, i ) != 0 )
		return 1;

	return 0;
    }

    public static int getRank( IntegerMatrix matrix ) {

	return GuIntegerMatrix.rankModuloTwo( matrix );
    }

    public static IntegerMatrix getH( IntegerMatrix ahi ) {

	int genus = ahi.getNumCols()/2;

	IntegerMatrix H = new IntegerMatrix ( genus );
	if( debug )
	    System.out.println(" genus = " + genus );
	for( int i=0; i<genus; i++)
	    for( int j=0; j< genus; j++)
		{
		    H.set( i, j, ahi.re[i][j+genus] );
		}
	return H;
    }

    
    private static IntegerVector getBasicEdges( CanonicalBasis canonicalBasis ) {
	IntegerVector basEdges = new IntegerVector ( canonicalBasis.numOfCycles );
	
	for(int j=0, k=0; j<canonicalBasis.numOfCycles; j++) {	    
	    k++;	    
	    while( canonicalBasis.statusOfEdges.re[k] != 0 ) 
		k++;
	    basEdges.set( j, k );
	}

	return basEdges;
    }

    static IntegerMatrix getTau( final CanonicalBasis canonicalBasis, int nl, int nr, int nb, IntegerVector tau0 ) {

	IntegerMatrix oldT = new IntegerMatrix ( canonicalBasis.T );

	IntegerVector edgeBranch = new IntegerVector ( canonicalBasis.edgeBranchPoint );
	IntegerVector edgeEnd    = new IntegerVector ( canonicalBasis.   edgeEndPoint );
	IntegerVector edgeStart  = new IntegerVector ( canonicalBasis. edgeStartPoint );

	/* We calculte the basic edges */

	IntegerVector basEdges = getBasicEdges( canonicalBasis );

	/* We move the origin from lambda+i0 to lambda-i0 */

	for ( int i = 1; i < nl ; i++)
	    for ( int j = 0; j < i ; j++)		
		canonicalBasis.branchPointPassedCutOfBranchPointClockWise( i, j );

	for ( int i = nl + nb + nr -2 ; i >= nl + nb ; i-- )
	    for ( int j =  nl + nb + nr - 1  ; j > i ; j-- )
		canonicalBasis.branchPointPassedCutOfBranchPointCounterClockWise( i, j );

	/* We reset the matrix T to the unit matrix */

	canonicalBasis.T.assignId();

	/* We find the corresponding edges after antiinvolution */	

	IntegerVector edgeCorrespondence = new IntegerVector ( canonicalBasis.numOf_g_Edges );
	
	for( int i=0; i<canonicalBasis.numOf_g_Edges; i++ ) {
	    int k        = edgeBranch.get(i);  
	    int beta     = edgeEnd.get(i);
	    int alpha    = edgeStart.get(i);
	    int tauBeta  = tau0.re[beta];
	    int tauAlpha = tau0.re[alpha];

	    int tauK     = k < nl || ( k >= nl + nb && k < nl + nb + nr ) ?
		k :
		canonicalBasis.numOfBranchPoints + nl - 1 - k;

	    edgeCorrespondence.set(i, canonicalBasis.gl.re[tauK][tauBeta] );

	    if( edgeCorrespondence.get(i) < 0 || canonicalBasis.g.re[tauK][tauBeta] != tauAlpha )
		throw new IllegalArgumentException(" The curve is not real " );
	}

	/* We move the origin back from lambda-i0 to lambda+i0 */

	for( int i = nl-1; i > 0 ; i-- )
	    for( int j = i - 1 ; j >= 0 ; j-- )
		canonicalBasis.branchPointPassedCutOfBranchPointCounterClockWise( i, j );

	for( int i = nl + nb ; i < nl + nb + nr -1   ; i++ )
	    for ( int j =  i + 1  ; j < nl + nb + nr ; j++ )
		canonicalBasis.branchPointPassedCutOfBranchPointClockWise( i, j );

	/* We calculate the antiinvolution on edges */

	IntegerMatrix antiinvOnEdges = new IntegerMatrix ( canonicalBasis.numOf_g_Edges );

	for( int i=0; i<canonicalBasis.numOf_g_Edges; i++ ) {
	    for( int j=0; j<canonicalBasis.numOf_g_Edges; j++ ) {
		int k = edgeCorrespondence.get(i);

		antiinvOnEdges.set( i, j, -canonicalBasis.T.re[k][j] );
	    }
	}

	/* We express the result of intiiinvolution on edges in terms of 
	   independent edges using 2-cell boundaries 
	*/

	IntegerMatrix antiinvOnIndEdges = new IntegerMatrix ( canonicalBasis.numOf_g_Edges );

	for( int i=0; i<canonicalBasis.numOf_g_Edges; i++ ) {
	    for( int j=0; j<canonicalBasis.numOf_g_Edges; j++ ) {
		int sum=0;
		for( int k=0; k < canonicalBasis.numOf_g_Edges; k++ )
		    sum += antiinvOnEdges.get(i,k) * canonicalBasis.dependenceMatrix.re[k][j];
		antiinvOnIndEdges.set( i, j, sum );
	    }
	}

	//	canonicalBasis.dependenceMatrix.print("Dep matr ");

	/* We recover the original value of T */

	canonicalBasis.T.assign( oldT );

	/* We extract the transformation to canonical basis from the 
	   vectors initialACycles, initialBCycles
	   We simultaneously calculate the action of antiinvolution in 
	   the original basis 
	*/

	IntegerMatrix toSimplectic = new IntegerMatrix ( canonicalBasis.numOfCycles );
	IntegerMatrix newTau       = new IntegerMatrix ( canonicalBasis.numOfCycles );

	for( int i=0; i<canonicalBasis.numOfCycles; i++) {
	    for( int j=0; j<canonicalBasis.numOfCycles; j++) {
		int k = basEdges.re[j];

		if(i < canonicalBasis.genus) 
		    toSimplectic.set( i, j, canonicalBasis.initialACycles.re[i                     ][k] );
		else 
		    toSimplectic.set( i, j, canonicalBasis.initialBCycles.re[i-canonicalBasis.genus][k] );

		int sum=0;

		for ( int k1=0; k1 < canonicalBasis.numOf_g_Edges; k1++)
		    sum += canonicalBasis.firstCycleBasis.re[i][k1] * antiinvOnIndEdges.get( k1, k );

		newTau.set( i, j, sum );
	    }
	}

	// newTau.print(" Antiinvolution in original basis " );

	/* we calculate the antiinvolution in simplectic basis */

	IntegerMatrix toSimplecticInv = GuIntegerMatrix.invert( toSimplectic );
	antiinvOnEdges.assignTimes( newTau, toSimplecticInv );
	newTau.assignTimes( toSimplectic, antiinvOnEdges );
	newTau.assignTranspose();

	// toSimplecticInv.print(" To simplectic Inverse");

	if( debug )
	    newTau.print("tau in simplectic basis");

	return newTau;
    }


    public static IntegerMatrix getAntiholomorphicInvolution
	( final CanonicalBasis canonicalBasis, int nl, int nr, int nb, IntegerVector tau0 ) {

	IntegerMatrix tau = getTau( canonicalBasis, nl, nr, nb, tau0 );

	return getAntiholomorphicInvolution( tau, getTransformToSymmetricBasis( tau ) );
    }


    static IntegerMatrix getSymmetricACycles( CanonicalBasis canonicalBasis, IntegerMatrix transf ) {

	/* Now we calculate the symmtrized basis */
	IntegerMatrix symmetricACycles = new IntegerMatrix ( canonicalBasis.genus , canonicalBasis.numOfEdges);
	
	for( int i=0; i < canonicalBasis.genus; i++)  {
	    for( int j=0; j < canonicalBasis.numOfEdges; j++) {
		
		int sum=0;

		for( int k=0; k < canonicalBasis.numOfCycles; k++) {
		    int k1 = k < canonicalBasis.genus ? 
			canonicalBasis.initialACycles.re[k                     ][j]:
			canonicalBasis.initialBCycles.re[k-canonicalBasis.genus][j];

		    sum += transf.get(k,i) * k1;
		}

		symmetricACycles.re[i][j] = sum;
	    }
	}

	return symmetricACycles;
    }


    static IntegerMatrix getSymmetricBCycles( CanonicalBasis canonicalBasis, IntegerMatrix transf ) {

	/* Now we calculate the symmtrized basis */
	IntegerMatrix symmetricBCycles = new IntegerMatrix ( canonicalBasis.genus , canonicalBasis.numOfEdges);
	
	for( int i=0; i<canonicalBasis.genus; i++)  {
	    for( int j=0; j< canonicalBasis.numOfEdges; j++) {
		
		int sum=0;
		
		for( int k = 0; k < canonicalBasis.numOfCycles; k++) {
		    int k1 = k < canonicalBasis.genus ?
			canonicalBasis.initialACycles.re[k                     ][j]:
			canonicalBasis.initialBCycles.re[k-canonicalBasis.genus][j];

		    sum += transf.get(k,i+canonicalBasis.genus) * k1;
		}

		symmetricBCycles.re[i][j]=sum;

	    }
	}

	return symmetricBCycles;
    }

    public static int getRealLocus( CanonicalBasis canonicalBasis, int nl, int nr, int nb, IntegerVector tau0 ) {

	/* we check if the real locus is empty */
  
	int realLocus=0;

	IntegerVector monodr = new IntegerVector ( canonicalBasis.numOfSheets );

	for( int k=0; k<canonicalBasis.numOfSheets; k++) {
	    monodr.re[k]=k;
	    if( monodr.re[k] == tau0.re[k] ) 
		realLocus=1;  
	}

	IntegerVector monodr2 = new IntegerVector ( canonicalBasis.numOfSheets );

	if( realLocus==0 ) {
	    for( int i = nl - 1 ; i >= 0; i--) {
		
		for( int k=0; k<canonicalBasis.numOfSheets; k++)
		    monodr2.re[k]=monodr.re[k];

		for( int k=0; k<canonicalBasis.numOfSheets; k++)
		    monodr.re[k]=monodr2.re[ canonicalBasis.g.re[i][k] ];

		if( debug ) {
		    System.out.println(" i = " + i );
		    // monodr.print( "  the current monodromy " );
		}

		for( int k=0; k<canonicalBasis.numOfSheets; k++)
		    if(monodr.re[k] == tau0.re[k] ) {
			realLocus=1; 
			// System.out.println(" The sheet #" + k + " is invariant.");
		    }
	    }
	}

  
	if( realLocus == 0 ) {
	    
	    for( int k=0; k<canonicalBasis.numOfSheets; k++)
		monodr.re[k]=k;

	    for( int i = nl + nb ; i < nl + nb + nr; i++) { 
		
		for( int k=0; k<canonicalBasis.numOfSheets; k++)
		    monodr2.re[k]=monodr.re[k];
		
		for( int k=0; k<canonicalBasis.numOfSheets; k++)
		    monodr.re[k] = canonicalBasis.g.re[i][ monodr2.re[k] ];
		
		if( debug ) {
		    System.out.println(" i = " + i );
		    // monodr.print( "  the current monodromy " );
		}

		for( int k=0; k<canonicalBasis.numOfSheets; k++)
		    if(monodr.re[k] == tau0.re[k] ) {
			realLocus=1; 
			// System.out.println(" The sheet #" + k + " is invariant.");
		    }
	    }
	}
    
	if( debug ) {

	    System.out.println(" realLocus = " + realLocus );

	    // toSimplectic.print("  G " );
	    // toSimplecticInv.print("  G^-1 " );
	}
	
	return realLocus;
    }

    
    public SymmetricBasis( SymmetricBasis aSymmetricBasis ) {
	this.assign( aSymmetricBasis );
    }

    public SymmetricBasis( IntegerMatrix monodrom, IntegerMatrix singular, IntegerVector disting,
			   int NL, int NR, int NB, IntegerVector aTau0 ) {

	super( monodrom, singular, disting );
	
	tau0 = new IntegerVector ( aTau0 );

	nl = NL;
	nr = NR;
	nb = NB;

	IntegerMatrix tau    = getTau( this, nl, nr, nb, tau0 );

	// tau.times( tau ).print( " tau^2 " ); 

	IntegerMatrix transf = getTransformToSymmetricBasis( tau );
	
	IntegerMatrix ahi    = getAntiholomorphicInvolution( tau, transf );

	IntegerMatrix symmetricACycles = getSymmetricACycles( this, transf );
	IntegerMatrix symmetricBCycles = getSymmetricBCycles( this, transf );

	H = getH( ahi );
	
	realLocus = getRealLocus( this, nl, nr, nb, tau0 );
	
	if( debug ) {
	    initialACycles.print("unsymmetric initial a-cycles");
	    initialBCycles.print("unsymmetric initial b-cycles");

	    symmetricACycles.print("symmetric initial a-cycles");
	    symmetricBCycles.print("symmetric initial b-cycles");
	}

	initialACycles = symmetricACycles;
	initialBCycles = symmetricBCycles;
    }

    static void test( CanonicalBasis canonicalBasis, int nl, int nr, int nb, IntegerVector tau0 ) {

	boolean  saveDebug = debug;

	debug = true;

	IntegerMatrix tau    = getTau( canonicalBasis, nl, nr, nb, tau0 );

	// tau.times( tau ).print( " tau^2 " ); 

	IntegerMatrix transf = getTransformToSymmetricBasis( tau );

	IntegerMatrix ahi    = getAntiholomorphicInvolution( tau, transf );

	IntegerMatrix symmetricACycles = getSymmetricACycles( canonicalBasis, transf );
	IntegerMatrix symmetricBCycles = getSymmetricBCycles( canonicalBasis, transf );

	ahi.print("Anitholomorphic involution in the Symmetric basis");
	symmetricACycles.print(" New A - cycles ");
	symmetricBCycles.print(" New B - cycles ");

	/*  test for simplecticity of new basis  */
	
	IntegerVector basEdges = getBasicEdges( canonicalBasis );
	
	IntegerMatrix antiinvOnEdges = new IntegerMatrix ( canonicalBasis.numOfCycles, canonicalBasis.numOfCycles );
	
	for( int i=0; i < canonicalBasis.numOfCycles; i++) {
	    for ( int j=0; j < canonicalBasis.numOfCycles; j++) {

		int sum=0;
      
		for ( int k=0; k < canonicalBasis.numOfCycles; k++) {
		    for( int k1=0; k1 < canonicalBasis.numOfCycles; k1++) {

			int alpha = basEdges.re[k];
			int beta  = basEdges.re[k1];

			int tauAlpha = i < canonicalBasis.genus ? 
			    symmetricACycles.re[i                     ][alpha] :
			    symmetricBCycles.re[i-canonicalBasis.genus][alpha];

			int tauBeta = j < canonicalBasis.genus ?
			    symmetricACycles.re[j                     ][beta] :
			    symmetricBCycles.re[j-canonicalBasis.genus][beta];

			sum += tauAlpha*tauBeta*canonicalBasis.firstCycleIntersections.re[k][k1];
		    }
		}
		
		antiinvOnEdges.set( i, j, sum );
	    }
	}

	antiinvOnEdges.print( " Intersection form in new basis " );

	int realLocus = getRealLocus( canonicalBasis, nl, nr, nb, tau0 );

	IntegerMatrix H = getH( ahi );

	H.print( "H" );

	int rankOfH    = GuIntegerMatrix.rankModuloTwo( H );
	int aInvariant = getInvariantA ( H );

	System.out.println( "\n  aInavriant = " + aInvariant + "%i  rank = " + rankOfH );

	debug = saveDebug;
	    
    }


}
