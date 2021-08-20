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

import java.io.Serializable;
import java.util.ArrayList;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.blas.IntegerMatrix;
import de.jtem.blas.IntegerVector;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.surface.DistinguishedPoint;
import de.jtem.riemann.surface.PointPassedCutEvent;
import de.jtem.riemann.surface.SingularPoint;

/** 
    These variables are defined during instantiation and stay unchanged afterwards:</br>

    <i>numOfSheets</i><br>
    <i>numOfBranchPoints</i><br>
    <i>numOfSingularPoints</i><br>
    <i>numOfEdges</i><br>
    <i>numOf_g_Edges</i><br>
    <i>numOf_h_Edges</i><br>
    <i>numOfPaths</i><br>
    <i>numOfDistinguishedPoints</i><br>
    <i>numOfCells</i><br>
    <br>
    <br>
    These variables are used only during the instantiation:<br> 
    <i>monodromyAboutInfinity</i><br>
    <br>
    <br>

    These variables are recalculated after each transformation, which are movements of
    the branch points:<br>
    
    <i>IntegerMatrix g  </i><br>
    <i>IntegerMatrix gl </i><br>
    <i>IntegerMatrix h  </i><br>
    <i>IntegerMatrix T  </i><br>
    <i>distinguishedPoints </i><br>
    <i>edgeStartPoint  </i><br>
    <i>edgeBranchPoint </i><br>
    <i>edgeEndPoint    </i><br>

  */
public class Transform implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    int numOfSheets;
    int numOfBranchPoints;
    int numOfSingularPoints;
    int numOfEdges;
    int numOf_g_Edges;
    int numOf_h_Edges;
    int numOfPaths;
    int numOfDistinguishedPoints;
    int numOfCells;

    final IntegerMatrix g  = new IntegerMatrix ();
    final IntegerMatrix gl = new IntegerMatrix ();
    final IntegerMatrix T  = new IntegerMatrix ();
    IntegerMatrix h;

    
    IntegerVector distinguishedPoints;

    final IntegerVector monodromyAboutInfinity = new IntegerVector ();
    final IntegerVector edgeStartPoint         = new IntegerVector ();
    final IntegerVector edgeBranchPoint        = new IntegerVector ();
    final IntegerVector edgeEndPoint           = new IntegerVector ();

    boolean debug = false;
    boolean passedCutDebug = false;

    private final IntegerVector ginv = new IntegerVector ();
    private final IntegerVector gtmp = new IntegerVector ();
    private final IntegerVector htmp = new IntegerVector ();


    public TransformState getState() {
	return new TransformState( this );
    }

    public void getState( TransformState aState ) {
	aState.save( this );
    }

    public void setState( TransformState aState ) {
	aState.load( this );
    }

    


    public void assign( Transform aTransform ) {

	numOfSheets              = aTransform.numOfSheets;
	numOfBranchPoints        = aTransform.numOfBranchPoints;
	numOfSingularPoints      = aTransform.numOfSingularPoints;
	numOfEdges               = aTransform.numOfEdges;
	numOf_g_Edges            = aTransform.numOf_g_Edges;
	numOf_h_Edges            = aTransform.numOf_h_Edges;
	numOfPaths               = aTransform.numOfPaths;
	numOfDistinguishedPoints = aTransform.numOfDistinguishedPoints;
	numOfCells               = aTransform.numOfCells;

	g .assign(  aTransform.g  );
	gl.assign(  aTransform.gl );
	T .assign(  aTransform.T  );

	if( aTransform.h == null ) h = null;
	else if( h == null ) h = new IntegerMatrix ( aTransform.h  );
	else h.assign(  aTransform.h  );

	if( aTransform.distinguishedPoints == null ) distinguishedPoints = null;
	else if( distinguishedPoints == null ) 
	    distinguishedPoints = new IntegerVector ( aTransform.distinguishedPoints  );
	else distinguishedPoints.assign(  aTransform.distinguishedPoints  );

	monodromyAboutInfinity.assign(  aTransform.monodromyAboutInfinity  );
	edgeStartPoint.assign(  aTransform.edgeStartPoint  );
	edgeBranchPoint.assign(  aTransform.edgeBranchPoint  );
	edgeEndPoint.assign(  aTransform.edgeEndPoint  );

	debug = aTransform.debug;
    }


    public Object clone() {
	return new Transform( this );
    }

    
    /**
       * Get the value of debug.
       * @return Value of debug.
       */
    public boolean getDebug() {return debug;}
    
    /**
       * Set the value of debug.
       * @param v  Value to assign to debug.
       */
    public void setDebug(boolean  v) {this.debug = v;}
    

    /**
       * Get the value of debug.
       * @return Value of debug.
       */
    public boolean getPassedCutDebug() {return passedCutDebug;}
    
    /**
       * Set the value of debug.
       * @param v  Value to assign to debug.
       */
    public void setPassedCutDebug(boolean  v) {this.passedCutDebug = v;}
    
    /**
       * Get the value of numOfSheets.
       * @return Value of numOfSheets.
       */
    public int getNumOfSheets() {return numOfSheets;}
    
    /**
       * Get the value of numOfBranchPoints.
       * @return Value of numOfBranchPoints.
       */
    public int getNumOfBranchPoints() {return numOfBranchPoints;}
    
    /**
       * Get the value of numOfSingularPoints.
       * @return Value of numOfSingularPoints.
       */
    public int getNumOfSingularPoints() {return numOfSingularPoints;}
    
     /**
       * Get the value of numOfEdges.
       * @return Value of numOfEdges.
       */
    public int getNumOfEdges() {return numOfEdges;}
    
    /**
       * Get the value of numOf_g_Edges.
       * @return Value of numOf_g_Edges.
       */
    public int getNumOf_g_Edges() {return numOf_g_Edges;}
    
    /**
     * Get the value of numOf_h_Edges.
     * @return Value of numOf_h_Edges.
     */
    public int getNumOf_h_Edges() {return numOf_h_Edges;}
    
    /**
       * Get the value of numOfPaths.
       * @return Value of numOfPaths.
       */
    public int getNumOfPaths() {return numOfPaths;}
    
    /**
       * Get the value of numOfDistinguishedPoints.
       * @return Value of numOfDistinguishedPoints.
       */
    public int getNumOfDistinguishedPoints() {return numOfDistinguishedPoints;}
    
    /**
       * Get the value of numOfCells.
       * @return Value of numOfCells.
       */
    public int getNumOfCells() {return numOfCells;}
    
    /**
       * Get the value of g.
       * @return Value of g.
       */
    public IntegerMatrix getG() {return g;}
   
    /**
       * Get the value of gl.
       * @return Value of gl.
       */
    public IntegerMatrix getGl() {return gl;}
    
    /**
       * Get the value of h.
       * @return Value of h.
       */
    public IntegerMatrix getH() {return h;}
    
    /**
       * Get the value of T.
       * @return Value of T.
       */
    public IntegerMatrix getT() {return T;}
       
    /**
       * Get the value of distinguishedPoints.
       * @return Value of distinguishedPoints.
       */
    public IntegerVector getDistinguishedPoints() {return distinguishedPoints;}

    /**
       * Get the value of monodromyAboutInfinity.
       * @return Value of monodromyAboutInfinity.
       */
    public IntegerVector getMonodromyAboutInfinity() {return monodromyAboutInfinity;}
    
    /**
       * Get the value of edgeStartPoint.
       * @return Value of edgeStartPoint.
       */
    public IntegerVector getEdgeStartPoint() {return edgeStartPoint;}
    
    /**
       * Get the value of edgeBranchPoint.
       * @return Value of edgeBranchPoint.
       */
    public IntegerVector getEdgeBranchPoint() {return edgeBranchPoint;}

    /**
       * Get the value of edgeEndPoint.
       * @return Value of edgeEndPoint.
       */
    public IntegerVector getEdgeEndPoint() {return edgeEndPoint;}

    protected Transform() {
    }
 
    public Transform( Transform aTransform ) {
	assign( aTransform );
    } 
    public Transform( IntegerMatrix monodrom ) {
	this( monodrom, null, null );
    } 
    public Transform( IntegerMatrix monodrom, IntegerMatrix singular ) {
	this( monodrom, singular, null );
    }  
    public Transform( IntegerMatrix monodrom, IntegerVector disting ) {
	this( monodrom, null, disting );
    }

    /**
     * Compute monodromy about infinity.
     * Methods assumes a COUNTERCLOCKWISE order of the branch points at this very moment.
     */
    private void computeMonodromyAboutInfinity( IntegerMatrix m ) {

	monodromyAboutInfinity.newSize( numOfSheets );
	
	for( int j=0; j<numOfSheets; j++ ) {
	    monodromyAboutInfinity.re[j] = j; 
	}

	for( int j=0; j<numOfSheets; j++ ) {
	    int current_g=j;
	    for( int i=0; i<numOfBranchPoints; i++ ) current_g=(int)m.re[i][current_g];
	    monodromyAboutInfinity.re[j]=current_g; 
	}

	if( debug )
	    monodromyAboutInfinity.print( "Monodromy about infinity" );
    }

    private void computeNumOfCells() {

	gtmp.assignZero();

	/* Now we calculte the number of infinte points == number of cells conaining infinite points */
	numOfCells=0;

	for( int j=0; j<numOfSheets; j++ ) {  
	    
	    if(gtmp.re[j]==0) {

		numOfCells++;
		int  current_g=j;	
		while( true ) {		
		    gtmp.re[current_g]=1;  
		    current_g=(int)monodromyAboutInfinity.re[current_g];
		    if(current_g==j) break;
		}
	    }
	}	
	
	if( debug )
	    System.out.println( "Number of Cells = " + numOfCells );
    }

    private void computeSystemOfEdges( IntegerMatrix m, IntegerMatrix s) {

	numOf_g_Edges = numOf_h_Edges = 0;

	g .newSize( numOfBranchPoints, numOfSheets );
	gl.newSize( numOfBranchPoints, numOfSheets );

	for( int i=0; i<numOfBranchPoints; i++ )
	    for( int j=0; j<numOfSheets; j++ ) {
		g.re[i][j]  = m.re[i][j]; 
		gl.re[i][j] = -2; 
	    }

	int current_edge=0;

	for( int i=0; i<numOfBranchPoints; i++ ) 
	    for( int j=0; j<numOfSheets; j++ ) {
		int current_g=(int)g.re[i][j];

		if(current_g==j) {
		    gl.re[i][j] = -1;

		} else if(gl.re[i][j]==-2) {
		    numOfCells++;
		    gl.re[i][j]=current_edge++; 
		    while(current_g!=j) {
			gl.re[i][current_g]=current_edge++;
			current_g=(int)g.re[i][current_g];
		    }
		}
	    }
    
    numOf_g_Edges = current_edge;

    if( s!=null ) {
	
	if( s.getNumCols() != numOfSheets )
	    throw new IllegalArgumentException( "matrices to not match" );
	    
	numOfSingularPoints = s.getNumRows();

	if( h == null )
	    h = new IntegerMatrix ( numOfSingularPoints, numOfSheets );
	else 
	    h.newSize( numOfSingularPoints, numOfSheets ); 

	for( int i=0; i<numOfSingularPoints; i++ ) 
	    for( int j=0; j<numOfSheets; j++ ) {			    
		if( s.re[i][j] == 0 )
		    h.re[i][j]=-1;
		else {
		    h.re[i][j]=current_edge++;
		    numOf_h_Edges++;
		}
	    }

	numOfEdges=current_edge;

	} else numOfSingularPoints = 0;

	numOfEdges = current_edge;
    }

    private void computeEdgeArrays() {

	edgeStartPoint .newSize( numOf_g_Edges );
	edgeBranchPoint.newSize( numOf_g_Edges );
	edgeEndPoint   .newSize( numOf_g_Edges );

	for( int i=0; i<numOfBranchPoints; i++ ) 
	    for( int j=0; j<numOfSheets; j++ ) {
		int current_g=(int)gl.re[i][j];
		if(current_g>-1) 
		    {
			int k=(int)g.re[i][j];
			edgeStartPoint.re [current_g]=j;     
			edgeBranchPoint.re[current_g]=i;  
			edgeEndPoint.re   [current_g]=k;  
		    }
	    }

	if( false || debug ) {
	    System.out.println("\n All g- Edges \n" );

	    edgeBranchPoint.print ( "Branch number = " );
	    edgeStartPoint.print( "Start point   = " );
	    edgeEndPoint.print( "End Point     = " );
	    
	    RuntimeException e = new RuntimeException( "stop ");

	    e.printStackTrace();

	    throw e;
	}
    }

    public Transform( IntegerMatrix monodrom, IntegerMatrix singular, IntegerVector disting ) {

	if( monodrom == null )
	    throw new IllegalArgumentException( "need an monodromie matrix" );

	numOfSheets       = monodrom.getNumCols();
	numOfBranchPoints = monodrom.getNumRows();

	ginv.newSize( numOfSheets );
	gtmp.newSize( numOfSheets );
	htmp.newSize( numOfSheets );

	computeMonodromyAboutInfinity( monodrom );
	computeNumOfCells();
	computeSystemOfEdges( monodrom, singular );
	computeEdgeArrays();

	if( disting != null ) {

	    numOfDistinguishedPoints = disting.size();

	    for( int i=0; i<numOfDistinguishedPoints; i++ ) 
		if( disting.re[i] < 0 || disting.re[i] >= numOfSheets )
		    throw new IllegalArgumentException( "distinguished points are out of range" );

	    distinguishedPoints = disting;

	} else numOfDistinguishedPoints = 0;

	numOfPaths = numOfEdges + numOfDistinguishedPoints;
	
	T.assignId( numOfPaths );
    }

    protected final void checkForBranchPointIndex( int anIndex ) {
	if( anIndex < 0 || anIndex >= numOfBranchPoints )
	    throw new IllegalArgumentException("branch point index "+anIndex+" is not in [0,"+numOfBranchPoints+"[");
    }

    protected final void checkForSingularPointIndex( int anIndex ) {
	if( anIndex < 0 || anIndex >= numOfSingularPoints )
	    throw new IllegalArgumentException("branch point index "+anIndex+" is not in [0,"+numOfSingularPoints+"[");
    }

    protected final void checkForDistinguishedPointIndex( int anIndex ) {
	if( anIndex < 0 || anIndex >= numOfDistinguishedPoints )
	    throw new IllegalArgumentException("branch point index "+anIndex+" is not in [0,"+numOfDistinguishedPoints+"[");
    }

    protected final void checkForBifurcation( int anIndex, int otherIndex ) {
	if( anIndex == otherIndex )
	    throw new IllegalArgumentException("not a bifurcation " + anIndex + " = " + otherIndex );
    }

    public void branchPointPassedCutOfBranchPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForBranchPointIndex( i_pnt );
	checkForBranchPointIndex( j_pnt );
	checkForBifurcation( i_pnt, j_pnt );

	if( passedCutDebug ) System.out.println("The Branch Point "+ i_pnt +" moved counter clockwise through Branch Point "+ j_pnt );
  
	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    int alpha=g.re[i_pnt][gamma]; 
	    ginv.re[alpha]=gamma;
	}

	/* The new monodromy matrix is calculated */
	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    int i1 = g.re[i_pnt][gamma]; 
	    i1 = g.re[j_pnt][i1]; 
	    gtmp.re[gamma]=ginv.re[i1] ;
	}

	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    g.re[j_pnt][gamma]=gtmp.re[gamma];
	}

	/* The new edge matrix is calculated */
	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    int alpha = g.re[i_pnt][gamma];
	    gtmp.re[gamma]=gl.re[j_pnt][alpha];
	}

	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    gl.re[j_pnt][gamma] = gtmp.re[gamma];
	}   

	/* The new matrix T is calculated */

	IntegerMatrix T1 = IntegerMatrix.id(  numOfPaths );

	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    int i1 = gl.re[j_pnt][gamma];

	    if(i1>-1) {
		int i2   = g.re[j_pnt][gamma];
		int beta = gl.re[i_pnt][i2];

		if(beta>-1)
		    T1.re[i1][beta]+=1;

		beta=gl.re[i_pnt][gamma];

		if(beta>-1)
		    T1.re[i1][beta]-=1;
	    }
	}
	
	T.assignTimes( T1 );

	computeEdgeArrays();
    }

    public int getEdgeNumber( int branchPointIndex, int sheetIndex ) {

	for( int i=0; i<numOf_g_Edges; i++ )

	    if( edgeBranchPoint.re[i] == branchPointIndex &&
		edgeStartPoint.re[i] == sheetIndex            )
		return i;

	throw new IllegalArgumentException();
    }

    public void branchPointPassedCutOfBranchPointClockWise( int i_pnt, int j_pnt ) {

	checkForBranchPointIndex( i_pnt );
	checkForBranchPointIndex( j_pnt );
	checkForBifurcation( i_pnt, j_pnt );

	if( passedCutDebug ) System.out.println("The Branch Point "+ i_pnt +" moved clockwise through Branch Point "+ j_pnt );

	/* The invesrse transfrom to g[i][*] is calculated */
  
	for( int alpha=0; alpha<numOfSheets; alpha++ )  {
	    int gamma = g.re[i_pnt][alpha]; ginv.re[gamma]=alpha;
	}

	/* The new monodromy matrix is calculated */
	for( int alpha=0; alpha<numOfSheets; alpha++ ) {
	    int gamma=g.re[i_pnt][alpha]; 
	    int i1=g.re[j_pnt][alpha]; 
	    gtmp.re[gamma]=g.re[i_pnt][i1] ;
	}

	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    g.re[j_pnt][gamma]=gtmp.re[gamma];
	}

	/* The new edge matrix is calculated */
	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    int alpha=ginv.re[gamma];
	    gtmp.re[gamma]=gl.re[j_pnt][alpha];
	}

	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    gl.re[j_pnt][gamma]=gtmp.re[gamma];
	}   

	/* The new matrix T is calculated */

	IntegerMatrix T1 = IntegerMatrix.id(  numOfPaths );
 
	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    int i1 = gl.re[j_pnt][gamma];

	    if(i1>-1) {
		int i2   = g.re[j_pnt][gamma];
		i2       = ginv.re[i2];
		int beta = gl.re[i_pnt][i2];

		if(beta>-1)
		    T1.re[i1][beta]-=1;

		i2   = ginv.re[gamma];
		beta = gl.re[i_pnt][i2];

		if(beta>-1)
		    T1.re[i1][beta]+=1;
	    }
	}

	T.assignTimes( T1 );

	computeEdgeArrays();

    }
    
    public void branchPointPassedCutOfSingularPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForBranchPointIndex( i_pnt );
	checkForSingularPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Branch Point "+ i_pnt +" moved counter clockwise through Singular Point "+ j_pnt );

	/* The new edge matrix is calculated */
	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
		int alpha      = g.re[i_pnt][gamma];
		htmp.re[gamma] = h.re[j_pnt][alpha];
	    }

	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    h.re[j_pnt][gamma] = htmp.re[gamma];
	}   
    }

    public void branchPointPassedCutOfSingularPointClockWise( int i_pnt, int j_pnt ) {

	checkForBranchPointIndex( i_pnt );
	checkForSingularPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Branch Point "+ i_pnt +" moved clockwise through Singular Point "+ j_pnt );

	/* The new edge matrix is calculated */
	for( int alpha=0; alpha<numOfSheets; alpha++ ) {
	    int gamma      = g.re[i_pnt][alpha];
	    htmp.re[gamma] = h.re[j_pnt][alpha];
	}

	for( int gamma=0; gamma<numOfSheets; gamma++ ) {
	    h.re[j_pnt][gamma] = htmp.re[gamma];
	}   
    }

    public void singularPointPassedCutOfBranchPointCounterClockWise( int i_pnt, int j_pnt ) {
	
	checkForSingularPointIndex( i_pnt );
	checkForBranchPointIndex( j_pnt );
	
	if( passedCutDebug ) System.out.println("The Singular Point "+ i_pnt +" moved counter clockwise through Branch Point "+ j_pnt );
	
	/* The new matrix T is calculated */

	IntegerMatrix T1 = IntegerMatrix.id( numOfPaths );

	for( int alpha=0; alpha<numOfSheets; alpha++ ) {
	    int i1 = gl.re[j_pnt][alpha];

	    if(i1>-1) {
		int i2   = g.re[j_pnt][alpha];
		int beta = h.re[i_pnt][i2];

		if(beta>-1) 
		    T1.re[i1][beta]+=1;

		beta = h.re[i_pnt][alpha];

		if(beta>-1) 
		    T1.re[i1][beta]-=1;
	    }
	}

	T.assignTimes( T1 );
    }

    public void singularPointPassedCutOfBranchPointClockWise( int i_pnt, int j_pnt ) {

	checkForSingularPointIndex( i_pnt );
	checkForBranchPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Singular Point "+ i_pnt +" moved clockwise through Branch Point "+ j_pnt );
	
	/* The new matrix T is calculated */
	
	IntegerMatrix T1 = IntegerMatrix.id(  numOfPaths );

	for( int alpha=0; alpha<numOfSheets; alpha++ ) {
	    int i1 = gl.re[j_pnt][alpha];

	    if(i1>-1) {
		int i2   = g.re[j_pnt][alpha];
		int beta = h.re[i_pnt][i2];

		if(beta>-1) 
		    T1.re[i1][beta]-=1;

		beta = h.re[i_pnt][alpha];

		if(beta>-1)
		    T1.re[i1][beta]+=1;
	    }
	}
	
	T.assignTimes( T1 );
    }

    public void singularPointPassedCutOfSingularPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForSingularPointIndex( i_pnt );
	checkForSingularPointIndex( j_pnt );
	checkForBifurcation( i_pnt, j_pnt );

	if( passedCutDebug ) System.out.println("The Singular Point "+ i_pnt +" moved counter clockwise through Singular Point "+ j_pnt );
    }

    public void singularPointPassedCutOfSingularPointClockWise( int i_pnt, int j_pnt ) {

	checkForSingularPointIndex( i_pnt );
	checkForSingularPointIndex( j_pnt );
	checkForBifurcation( i_pnt, j_pnt );

	if( passedCutDebug ) System.out.println("The Singular Point "+ i_pnt +" moved clockwise through Singular Point "+ j_pnt );
    }

    public void branchPointPassedCutOfDistinguishedPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForBranchPointIndex( i_pnt );
	checkForDistinguishedPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Branch Point "+ i_pnt +" moved counter clockwise through Distinguished Point "+ j_pnt );
 
	/* The element g^-1[i][beta] is calculated */
 
	int alpha = -1;
	int beta  = distinguishedPoints.re[j_pnt];

	for( int i1=0; i1 < numOfSheets; i1 ++ ) {
	    int i2 = g.re[i_pnt][i1];

	    if(i2 == beta) {
		alpha = i1;
		break;
	    }
	}

	distinguishedPoints.re[j_pnt] = alpha;

	/* The new matrix T is calculated */

	IntegerMatrix T1 = IntegerMatrix.id(  numOfPaths );

	int gamma = gl.re[i_pnt][alpha];

	if(gamma>-1) 
	    T1.re[j_pnt+numOfEdges][gamma]-=1;
	
	T.assignTimes( T1 );
    }

    public void branchPointPassedCutOfDistinguishedPointClockWise( int i_pnt, int j_pnt ) {

	checkForBranchPointIndex( i_pnt );
	checkForDistinguishedPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Branch Point "+ i_pnt +" moved clockwise through Distinguished Point "+ j_pnt );

	/* The element beta=g[i][alpha] is calculated */
 
	int alpha = distinguishedPoints.re[j_pnt];
	int beta  = g.re[i_pnt][alpha];

	distinguishedPoints.re[j_pnt] = beta;

	/* The new matrix T is calculated */

	IntegerMatrix T1 = IntegerMatrix.id(  numOfPaths );
	
	int gamma = gl.re[i_pnt][alpha];

	if(gamma>-1) 
	    T1.re[j_pnt+numOfEdges][gamma]+=1;

	T.assignTimes( T1 );
    }

    public void singularPointPassedCutOfDistinguishedPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForSingularPointIndex( i_pnt );
	checkForDistinguishedPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Singular Point "+ i_pnt +" moved counter clockwise through Distinguished Point "+ j_pnt );
 
	int alpha = distinguishedPoints.re[j_pnt];

	/* The new matrix T is calculated */
	
	IntegerMatrix T1 = IntegerMatrix.id(  numOfPaths );

	int gamma = h.re[i_pnt][alpha];

	if(gamma>-1) 
	    T1.re[j_pnt+numOfEdges][gamma]-=1;

	T.assignTimes( T1 );	
    }

    
    public void singularPointPassedCutOfDistinguishedPointClockWise( int i_pnt, int j_pnt ) {

	checkForSingularPointIndex( i_pnt );
	checkForDistinguishedPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Singular Point "+ i_pnt +" moved clockwise through Distinguished Point "+ j_pnt );
 
	int alpha = distinguishedPoints.re[j_pnt];

	/* The new matrix T is calculated */
	
	IntegerMatrix T1 = IntegerMatrix.id(  numOfPaths );

	int gamma=h.re[i_pnt][alpha];

	if(gamma>-1) 
	    T1.re[j_pnt+numOfEdges][gamma]+=1;

	T.assignTimes( T1 );	
    }

    public void distinguishedPointPassedCutOfBranchPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForDistinguishedPointIndex( i_pnt );
	checkForBranchPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Distinguished Point "+ i_pnt +" moved counter clockwise through Branch Point "+ j_pnt );
    }
    
    public void distinguishedPointPassedCutOfBranchPointClockWise( int i_pnt, int j_pnt ) {

	checkForDistinguishedPointIndex( i_pnt );
	checkForBranchPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Distinguished Point "+ i_pnt +" moved clockwise through Branch Point "+ j_pnt );
    }

    public void distinguishedPointPassedCutOfSingularPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForDistinguishedPointIndex( i_pnt );
	checkForSingularPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Distinguished Point "+ i_pnt +" moved counter clockwise through Singular Point "+ j_pnt );
    }
    
    public void distinguishedPointPassedCutOfSingularPointClockWise( int i_pnt, int j_pnt ) {

	checkForDistinguishedPointIndex( i_pnt );
	checkForSingularPointIndex( j_pnt );

	if( passedCutDebug ) System.out.println("The Distinguished Point "+ i_pnt +" moved clockwise through Singular Point "+ j_pnt );
    }

    public void distinguishedPointPassedCutOfDistinguishedPointCounterClockWise( int i_pnt, int j_pnt ) {

	checkForDistinguishedPointIndex( i_pnt );
	checkForDistinguishedPointIndex( j_pnt );
	checkForBifurcation( i_pnt, j_pnt );

	if( passedCutDebug ) System.out.println("The Distinguished Point "+ i_pnt +" moved counter clockwise through Distinguished Point "+ j_pnt );
    }
    
    public void distinguishedPointPassedCutOfDistinguishedPointClockWise( int i_pnt, int j_pnt ) {

	checkForDistinguishedPointIndex( i_pnt );
	checkForDistinguishedPointIndex( j_pnt );
	checkForBifurcation( i_pnt, j_pnt );

	if( passedCutDebug ) System.out.println("The Distinguished Point "+ i_pnt +" moved clockwise through Distinguished Point "+ j_pnt );
    }

    public void process( PointPassedCutEvent e ) { 

	int pointIndex = e.point.getIndex();
	int cutIndex   = e.cut.getIndex();

	if( e.clockWise ) {

	    if( e.point instanceof BranchPoint  ) {

		if( e.cut instanceof BranchPoint ) {
		    branchPointPassedCutOfBranchPointClockWise( pointIndex, cutIndex );
		} else if( e.cut instanceof SingularPoint ) {
		    branchPointPassedCutOfSingularPointClockWise( pointIndex, cutIndex );
		} else if( e.cut instanceof DistinguishedPoint ) {
		    branchPointPassedCutOfDistinguishedPointClockWise( pointIndex, cutIndex );
		}
	    } else if( e.point instanceof SingularPoint ) {
		
		if( e.cut instanceof BranchPoint ) {
		    singularPointPassedCutOfBranchPointClockWise( pointIndex, cutIndex );
		} else if( e.cut instanceof SingularPoint ) {
		    singularPointPassedCutOfSingularPointClockWise( pointIndex, cutIndex );

		} else if( e.cut instanceof DistinguishedPoint ) {
		    singularPointPassedCutOfDistinguishedPointClockWise( pointIndex, cutIndex );
		}
	    } else if( e.point instanceof DistinguishedPoint ) {
		    
		    if( e.cut instanceof BranchPoint ) {
			distinguishedPointPassedCutOfBranchPointClockWise( pointIndex, cutIndex );
		    } else if( e.cut instanceof SingularPoint ) {
			distinguishedPointPassedCutOfSingularPointClockWise( pointIndex, cutIndex );
		    } else if( e.cut instanceof DistinguishedPoint ) {
			distinguishedPointPassedCutOfDistinguishedPointClockWise( pointIndex, cutIndex );
		    }
		} 
	} else {

	    if( e.point instanceof BranchPoint  ) {

		if( e.cut instanceof BranchPoint ) {
		    branchPointPassedCutOfBranchPointCounterClockWise( pointIndex, cutIndex );
		} else if( e.cut instanceof SingularPoint ) {
		    branchPointPassedCutOfSingularPointCounterClockWise( pointIndex, cutIndex );
		} else if( e.cut instanceof DistinguishedPoint ) {
		    branchPointPassedCutOfDistinguishedPointCounterClockWise( pointIndex, cutIndex );
		}
	    } else if( e.point instanceof SingularPoint ) {
		
		if( e.cut instanceof BranchPoint ) {
		    singularPointPassedCutOfBranchPointCounterClockWise( pointIndex, cutIndex );
		} else if( e.cut instanceof SingularPoint ) {
		    singularPointPassedCutOfSingularPointCounterClockWise( pointIndex, cutIndex );

		} else if( e.cut instanceof DistinguishedPoint ) {
		    singularPointPassedCutOfDistinguishedPointCounterClockWise( pointIndex, cutIndex );
		}
	    } else if( e.point instanceof DistinguishedPoint ) {
		    
		    if( e.cut instanceof BranchPoint ) {
			distinguishedPointPassedCutOfBranchPointCounterClockWise( pointIndex, cutIndex );
		    } else if( e.cut instanceof SingularPoint ) {
			distinguishedPointPassedCutOfSingularPointCounterClockWise( pointIndex, cutIndex );
		    } else if( e.cut instanceof DistinguishedPoint ) {
			distinguishedPointPassedCutOfDistinguishedPointCounterClockWise( pointIndex, cutIndex );
		    }
	    }

	}
	    
    }

    public void process( ArrayList listOfPasses ) {

	int numOfPasses = listOfPasses.size();

	for( int i=0; i<numOfPasses; i++ ) {
	    process( (PointPassedCutEvent)listOfPasses.get(i) );
	}
    }
	

    public void print() {

	g.print( " Matrix g " );
	gl.print( " Matrix gl " );
	h.print( " Matrix h " );
	distinguishedPoints.print( " Array of Distinguished points = " );
	T.print( " Matrix T " );		    
    }

    public void computeEdgeIntegrals( ComplexMatrix integrals, ComplexVector edgeIntegrals ) {
	computeEdgeIntegrals( integrals, null, (ComplexVector)null, edgeIntegrals );
    }

    public void computeEdgeIntegrals( ComplexMatrix integrals, ComplexMatrix sing_integrals, ComplexVector edgeIntegrals ) {
	computeEdgeIntegrals( integrals, sing_integrals, (ComplexVector)null, edgeIntegrals );
    }

    public void computeEdgeIntegrals( ComplexMatrix integrals, ComplexMatrix sing_integrals, ComplexMatrix distinguishedIntegrals, ComplexVector edgeIntegrals ) {
	ComplexVector dInt = null;

	if( numOfDistinguishedPoints > 0 ) {
	    if(    distinguishedIntegrals.getNumRows() != numOfDistinguishedPoints
		   || distinguishedIntegrals.getNumCols() != numOfSheets ) 
		throw new IllegalArgumentException( "distinguishedIntegrals has wrong size" );

	    dInt = new ComplexVector( numOfDistinguishedPoints );
  
	    for(int i1=0; i1< numOfDistinguishedPoints; i1++) {
		int  i2 = distinguishedPoints.re[i1];

		dInt.re[i1] = distinguishedIntegrals.re[i1][i2];
		dInt.im[i1] = distinguishedIntegrals.im[i1][i2];
	    }
	}

	computeEdgeIntegrals( integrals, sing_integrals, dInt, edgeIntegrals );
    }
    
    public void computeEdgeIntegrals( ComplexMatrix integrals, ComplexMatrix sing_integrals, ComplexVector distinguishedIntegrals, ComplexVector edgeIntegrals ) {

	if( integrals == null ) 
	    throw new IllegalArgumentException( "no matrix" );

	if(    integrals.getNumRows() != numOfBranchPoints
	    || integrals.getNumCols() != numOfSheets ) 
	    throw new IllegalArgumentException( "integrals matrix has wrong size" );

	if( edgeIntegrals == null )
	    throw new IllegalArgumentException( "no data grid for the edgeIntegrals" );

	if( sing_integrals == null ) {

	    if( numOfSingularPoints != 0 )
		throw new IllegalArgumentException( "no matrix of integrals around singular points is defined but such points exist" );

	} else if(    sing_integrals.getNumRows() != numOfSingularPoints
		   || sing_integrals.getNumCols() != numOfSheets ) 
	    throw new IllegalArgumentException( "singular integrals matrix has wrong size" );

	if( distinguishedIntegrals == null ) {

	    if( numOfDistinguishedPoints != 0 ) 
		throw new IllegalArgumentException( "no vector of integrals to distinguished points is defined but such points exist" );

	} else if( distinguishedIntegrals.size() != numOfDistinguishedPoints )
	    throw new IllegalArgumentException( "distinguishedIntegrals array has wrong size" );

	edgeIntegrals.newSize( numOfPaths );
	
  
	for( int i1=0; i1<numOfPaths; i1++ ) {
	    double xr=0;
	    double xi=0; 

	    for( int i2=0; i2<numOfBranchPoints; i2++ )
		for( int i3=0; i3<numOfSheets; i3++ ) {		
		    int i4 = gl.re[i2][i3];
		    
		    if( i4 > -1 ) {
			xr += T.re[i1][i4] * integrals.re[i2][i3];  
			xi += T.re[i1][i4] * integrals.im[i2][i3];
		    }
		}   
       
	    if( sing_integrals != null )
 
		for( int i2=0; i2<numOfSingularPoints; i2++ )
		    for( int i3=0; i3<numOfSheets; i3++ ) {		
			int i4 = h.re[i2][i3];

			if( i4 >- 1 ) {
			    xr += T.re[i1][i4] * sing_integrals.re[i2][i3];
			    xi += T.re[i1][i4] * sing_integrals.im[i2][i3];
			}
		    }

	    if( distinguishedIntegrals != null )

		for( int i2=0; i2<numOfDistinguishedPoints; i2++ )  {
		    int i4 = i2+numOfEdges;

		    xr += T.re[i1][i4] * distinguishedIntegrals.re[i2];  
		    xi += T.re[i1][i4] * distinguishedIntegrals.im[i2];
		}


	    edgeIntegrals.re[i1] = xr; 
	    edgeIntegrals.im[i1] = xi;
	}
    }

}



