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

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import de.jtem.blas.ComplexVector;
import de.jtem.blas.IntegerVector;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.homologie.Transform;
import de.jtem.riemann.util.RadialPathGenerator;

public abstract class RamifiedCovering  implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    public static final double EPS = 1e-12;
    
    protected BranchPoint[]               branchPoint;
    protected SingularPoint[]           singularPoint;
    protected DistinguishedPoint[] distinguishedPoint;

    protected int numOfBranchPoints           = 0;
    protected int numOfSingularPoints         = 0;  
    protected int numOfDistinguishedPoints    = 0;
    protected int numOfSurfacePoints = 0;
    
    protected Origin origin;
   
    protected boolean uptodate      = false;    
    protected boolean isInitialized = false;

    private boolean enablePropertyChange      = true;
    private boolean propertyChangeWasNotFired = false;

    protected RadialPathGenerator pathGenerator;
    
    ArrayList listOfPasses = new ArrayList( numOfSurfacePoints );

    private Complex tmp = new Complex();

    protected PropertyChangeSupport 
	propertyChangeSupport = new PropertyChangeSupport( this );


	protected SurfacePoint[] initSurfacePoint;

	protected SurfacePoint[] surfacePoint;

    public void addPropertyChangeListener( PropertyChangeListener listener ) {
	
	propertyChangeSupport.addPropertyChangeListener( listener);
    }

    public void removePropertyChangeListener( PropertyChangeListener listener ) {
	propertyChangeSupport.removePropertyChangeListener( listener );
    }

    public boolean getEnablePropertyChange() {
	return enablePropertyChange;
    }

    public void setEnablePropertyChange( boolean aBool ) {

	if( aBool == enablePropertyChange )
	    return;

	enablePropertyChange = aBool;
	
	if( aBool ) {
	    if( propertyChangeWasNotFired )
		propertyChangeSupport.firePropertyChange( null, null, null );
	} else {
	    propertyChangeWasNotFired = false;
	}
    }

    public void firePropertyChange( String propertyName) {
	firePropertyChange( propertyName, enablePropertyChange );
    }

    public void firePropertyChange( String propertyName, boolean enablePropertyChange ) {
	if( enablePropertyChange ) 
	    propertyChangeSupport.firePropertyChange( propertyName, null, null );
	
	else
	    propertyChangeWasNotFired = true;
    }
    public Origin getOrigin() {return origin;}
    
    public abstract Transform getTransform();

    public BranchPoint[] getBranchPoints() {
	return branchPoint;
    }

    public SingularPoint [] getSingularPoints() {
	return singularPoint;
    }

    public DistinguishedPoint [] getDistinguishedPoints() {
	return distinguishedPoint;
    }

    public RamifiedCoveringState getState() {
	return new RamifiedCoveringState( this );
    }

    public void getState( RamifiedCoveringState aState ) {
	aState.save( this );
    }

    public void setState( RamifiedCoveringState aState ) {
	aState.load( this );

	outdate();

	firePropertyChange( "state");
    }

    public boolean getUptodate() {return uptodate;}

    public void setUptodate( boolean aBool) {
	if( aBool )
	    update();
	else
	    uptodate = false;
    }
    
    public void initCompute() {

	numOfBranchPoints        = branchPoint        == null ? 0 : branchPoint       .length;
	numOfSingularPoints      = singularPoint      == null ? 0 : singularPoint     .length;
	numOfDistinguishedPoints = distinguishedPoint == null ? 0 : distinguishedPoint.length;
	
	SurfacePoint [] pointsToAvoid = getPointsToAvoid();
	
	pathGenerator = new RadialPathGenerator( pointsToAvoid, origin );
	
	init( branchPoint );
	init( singularPoint );
	init( distinguishedPoint );

	initSurfacePoint = getSurfacePoints();
	surfacePoint     = getSurfacePoints();
	
	Arrays.sort( initSurfacePoint );
	Arrays.sort( surfacePoint );

	outdate();

	isInitialized = true;
    }

    public IntegerVector getDistinguishedVector() {

	if( numOfDistinguishedPoints < 1 )
	    return null;

	IntegerVector anIntegerVector =  new IntegerVector ( numOfDistinguishedPoints );

	for( int i=0; i<numOfDistinguishedPoints; i++ ) {
	    anIntegerVector.re[i] = distinguishedPoint[i].sheetNumber;
	}

	return anIntegerVector;
    }

    public void outdate() {
	uptodate = false;	
    }

    public void update() {

	if( !isInitialized )
	    initCompute();

	if( uptodate )
	    return;

	if( getTransform() != null )
	    getTransform().process( listOfPasses );

	listOfPasses.clear();

	uptodate = true;
    }

    public int getNumOfBranchPoints() {return numOfBranchPoints;}

    public int getNumOfSingularPoints() {return numOfSingularPoints;}

    public int getNumOfDistinguishedPoints() {return numOfDistinguishedPoints;}

    public BranchPoint getBranchPoint( int anIndex ) {
	return branchPoint[anIndex];
    }

    public SingularPoint getSingularPoint( int anIndex ) {
	return singularPoint[anIndex];
    }

    public DistinguishedPoint getDistinguishedPoint( int anIndex ) {
	return distinguishedPoint[anIndex];
    }

    public final BranchPoint getBranchPointWithCoords( Complex coords ) {
	return getBranchPointWithCoords( coords, EPS*EPS );
    }

    public final BranchPoint getBranchPointWithCoords( Complex coords, double epsSqr ) {

	for ( int i=0; i<branchPoint.length; i++ ) {
	    if( coords.distSqr( branchPoint[i] ) < epsSqr )
		return branchPoint[i];
	}
	
	return null;
    }
 
    public final int getIndexOfBranchPointWithCoords( Complex coords, double epsSqr ) {

	return SurfacePoint.getIndexOfSurfacePointWithCoords( branchPoint, coords, epsSqr );
    }

    public final int getIndexOfSingularPointWithCoords( Complex coords, double epsSqr ) {

	return SurfacePoint.getIndexOfSurfacePointWithCoords( singularPoint, coords, epsSqr );
    }

    public final int getIndexOfDistinguishedPointWithCoords( Complex coords, double epsSqr ) {

	return SurfacePoint.getIndexOfSurfacePointWithCoords( distinguishedPoint, coords, epsSqr );
    }

    public final int getIndexOfBranchPointWithCoords( Complex coords ) {

	return SurfacePoint.getIndexOfSurfacePointWithCoords( branchPoint, coords, EPS );
    }

    public final int getIndexOfSingularPointWithCoords( Complex coords ) {

	return SurfacePoint.getIndexOfSurfacePointWithCoords( singularPoint, coords, EPS );
    }

    public final int getIndexOfDistinguishedPointWithCoords( Complex coords ) {

	return SurfacePoint.getIndexOfSurfacePointWithCoords( distinguishedPoint, coords, EPS );
    }

    //public void remove( SurfacePoint thePoint, Field.Complex newCoords ) {
    //	thePoint.setNewCoords( newCoords );
    //} 

    protected SurfacePoint [] getSurfacePoints() {
	
	    numOfSurfacePoints = numOfBranchPoints + numOfSingularPoints + numOfDistinguishedPoints;
	
	SurfacePoint [] someSurfacePoint = new SurfacePoint[ numOfSurfacePoints ];
	
	    for( int i=0; i < numOfBranchPoints; i++ )
	        someSurfacePoint[ i ] = branchPoint[i];
	
	    for( int i=0; i < numOfSingularPoints; i++ )
	        someSurfacePoint[ i + numOfBranchPoints ] = singularPoint[i];
	
	    for( int i=0; i < numOfDistinguishedPoints; i++ )
	        someSurfacePoint[ i + numOfBranchPoints + numOfSingularPoints ] = distinguishedPoint[i];
	
	return someSurfacePoint;
	}

	protected SurfacePoint [] getPointsToAvoid() {
		
		int numOfPointsToAvoid = numOfBranchPoints + numOfSingularPoints;
		
		SurfacePoint [] someSurfacePoint = new SurfacePoint[ numOfPointsToAvoid ];
		
		for( int i=0; i < numOfBranchPoints; i++ )
			someSurfacePoint[ i ] = branchPoint[i];
		
		for( int i=0; i < numOfSingularPoints; i++ )
			someSurfacePoint[ i + numOfBranchPoints ] = singularPoint[i];	
		
		return someSurfacePoint;
	}
	
	public int getNumOfSurfacePoints() { return numOfSurfacePoints; }

	public final int getIndexOfSurfacePointWithCoords(Complex coords, double epsSqr) {
	
	return SurfacePoint.getIndexOfSurfacePointWithCoords( surfacePoint, coords, epsSqr );
	}

	public final int getIndexOfSurfacePointWithCoords(Complex coords) {
	
	return SurfacePoint.getIndexOfSurfacePointWithCoords( surfacePoint, coords, EPS );
	}

	public void computeLoops(ComplexVector[] branchPointLoop, ComplexVector[] singularPointLoop, ComplexVector[] distinguishedPointPath) {
	
		pathGenerator.update();
		
	pathGenerator.computeLoopsAroundTroublePoints( branchPoint, branchPointLoop );
	
	if( singularPoint != null )
		pathGenerator.computeLoopsAroundTroublePoints( singularPoint, singularPointLoop );
	
		//lambdaPlane.computeLoopsAroundPoints( singularPoint, singularPointLoop );
	
	if( distinguishedPoint != null )
		pathGenerator.computePathsToPoints( distinguishedPoint, distinguishedPointPath );
		//lambdaPlane.computePathsToPoints( distinguishedPoint, distinguishedPointPath );	
	}

	protected void init(SurfacePoint [] point) {
	if( point == null )
	    return;
	
	for( int i=0; i<point.length; i++ )
	    point[i].setSurface( this );
	
	Arrays.sort( point );	    
	
	for( int i=0; i<point.length; i++ )
		point[i].index = i;
	}
	
	

	private void switchSurfacePoints( int i, int j ) {

		SurfacePoint tmp = surfacePoint[i];
		
		surfacePoint[i] = surfacePoint[j];
		surfacePoint[j] = tmp;
	}
	
	public void move( SurfacePoint thePoint, Complex newCoords ) {

		double deltaArg = thePoint.getDeltaArgTo( newCoords );
		double newAbs   = thePoint.getAbsFor( newCoords );

		int indexOfThePoint      = SurfacePoint.getIndexOfSurfacePoint( surfacePoint, thePoint );
		int firstIndexOfThePoint = indexOfThePoint;

		if( thePoint.isToTheRightOf( newCoords ) ) { // move is counter clock wise
			
			while( true ) { 

				int indexOfOther = ( indexOfThePoint + 1 ) % numOfSurfacePoints;

				SurfacePoint other = surfacePoint[ indexOfOther ];

				if( !other.isToTheRightOf( newCoords ) || other.isToTheRightOf( thePoint ) )
					break;

				//System.out.println( "in loop: " + other.arg + " <=> " + thePoint.arg );
				
				double t = deltaArg == 0 ? 0.5 : thePoint.getDeltaArgTo( other ) / deltaArg;

				if( t < 0 || t > 1 )
					throw new RuntimeException("t = "+t+" is not valid");

				if( thePoint.abs * ( 1-t ) + newAbs * t < other.abs )
					listOfPasses.add( new PointPassedCutEvent( thePoint, other, 0, false ) );
				else
					listOfPasses.add( new PointPassedCutEvent( other, thePoint, 0, true  ) );


				switchSurfacePoints( indexOfThePoint, indexOfOther );

				indexOfThePoint = indexOfOther;

				if( indexOfThePoint == firstIndexOfThePoint )
					break;
			}

		} else {
			
			while( true ) { 

				int indexOfOther = ( indexOfThePoint + numOfSurfacePoints - 1 ) % numOfSurfacePoints;

				SurfacePoint other = surfacePoint[ indexOfOther ];

				if( other.isToTheRightOf( newCoords ) || !other.isToTheRightOf( thePoint ) )
					break;

				//System.out.println( "in loop: " + other.arg + " <=> " + thePoint.arg );

				double t = deltaArg == 0 ? 0.5 : thePoint.getDeltaArgTo( other ) / deltaArg;

				if( t < 0 || t > 1 )
					throw new RuntimeException("t = "+t+" is not valid");

				if( thePoint.abs * ( 1-t ) + newAbs * t < other.abs ) 
					listOfPasses.add( new PointPassedCutEvent( thePoint, other, 0, true ) );
				else
					listOfPasses.add( new PointPassedCutEvent( other, thePoint, 0, false  ) );

				switchSurfacePoints( indexOfThePoint, indexOfOther );

				indexOfThePoint = indexOfOther;

				if( indexOfThePoint == firstIndexOfThePoint )
					break;
			}
		}

		thePoint.setNewCoords(newCoords);
		
		outdate();
		firePropertyChange( "coords" );
	}
}













