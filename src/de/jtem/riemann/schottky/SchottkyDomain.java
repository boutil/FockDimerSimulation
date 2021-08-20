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

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.numericalMethods.geometry.meshGeneration.ruppert.Ruppert;

public class SchottkyDomain implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    
    Schottky schottky;
    
    Ruppert triangulator;

    double minAngle  = 10;
    double maxArea   = 1;
    
    double xDist = 1;
    double yDist = 1;

    double xMax, yMax, xMin, yMin;

    double boundary [][];

    double [] point;

    /**
     * Get the value of xMax.
     * @return Value of xMax.
     */
    public double getXMax() {
	return xMax;
    }

    /**
     * Get the value of yMax.
     * @return Value of yMax.
     */
    public double getYMax() {
	return yMax;
    }
      
    /**
     * Get the value of xMin.
     * @return Value of xMin.
     */
    public double getXMin() {
	return xMin;
    }

    /**
     * Get the value of yMin.
     * @return Value of yMin.
     */
    public double getYMin() {
	return yMin;
    }
    
    /**
     * Get the value of yDist.
     * @return Value of yDist.
     */
    public double getYDist() {
	return yDist;
    }
    
    /**
     * Set the value of yDist.
     * @param v  Value to assign to yDist.
     */
    public void setYDist(double  v) {
	if( this.yDist == v )
	    return;

	this.yDist = v;

	update();
    }
    
    /**
     * Get the value of xDist.
     * @return Value of xDist.
     */
    public double getXDist() {return xDist;}
    
    /**
     * Set the value of xDist.
     * @param v  Value to assign to xDist.
     */
    public void setXDist(double  v) {
	if( this.xDist == v )
	    return;
	
	this.xDist = v;

	update();
    }
    
    int circleDiscr = 12;
    
    /**
     * Get the value of circleDiscr.
     * @return Value of circleDiscr.
     */
    public int getCircleDiscr() {
	return circleDiscr;
    }
    
    /**
     * Set the value of circleDiscr.
     * @param v  Value to assign to circleDiscr.
     */
    public void setCircleDiscr(int  v) {
	if( this.circleDiscr == v )
	    return;
	
	this.circleDiscr = v;

	update();
    }
    
    /**
     * Get the value of maxArea.
     * @return Value of maxArea.
     */
    public double getMaxArea() {
	return maxArea;
    }
    
    /**
     * Set the value of maxArea.
     * @param v  Value to assign to maxArea.
     */
    public void setMaxArea(double  v) {

	final boolean rebuildTriangulator = v > this.maxArea;

	if( this.maxArea == v )
	    return;

	this.maxArea = v;
	
	updateTriangulator( rebuildTriangulator );
    }
    
    /**
     * Get the value of minAngle.
     * @return Value of minAngle.
     */
    public double getMinAngle() {
	return minAngle;
    }
    
    /**
     * Set the value of minAngle.
     * @param v  Value to assign to minAngle.
     */
    public void setMinAngle(double  v) {

	final boolean rebuildTriangulator = v < this.minAngle;

	if( this.minAngle == v )
	    return;
	
	this.minAngle = v;

	updateTriangulator( rebuildTriangulator );
    }
        
    /**
     * Get the value of schottky.
     * @return Value of schottky.
     */
    public Schottky getSchottky() {
	return schottky;
    }
    
    /**
     * Set the value of schottky.
     * @param v  Value to assign to schottky.
     */
    public void setSchottky(Schottky  v) {
	
	if( v == schottky )
	    return;
	
	this.schottky = v;
	
	update();
    }


    public SchottkyDomain() {	
    }
    
    public SchottkyDomain( final Schottky schottky ) {

	setSchottky( schottky );
    }

    public Ruppert getTriangulator() {
	return triangulator;
    }

    public static double [][] getBoundary( Schottky schottky, 
					   double xDist, double yDist, int circleDiscr ) {
	
	double [][] boundary = new double[ 2*schottky.getNumGenerators() + 1][];

	boundary[0] = getOuterBoundary( schottky, xDist, yDist );

	for( int i=0, k=1; i < schottky.getNumGenerators(); i++ ) {
	    
	    double [][] circlePair = getCirclePair( schottky, i, circleDiscr );
	    
	    boundary[k++] = circlePair[0];
	    boundary[k++] = circlePair[1];
	}

	return boundary;
    }

    public static double [] getOuterBoundary( Schottky schottky, double xDist, double yDist ) {

	double [] bound = SchottkyDomainSampler.getBound( schottky );

	double xMin = bound[0] - xDist + 1e-2;
	double yMin = bound[1] - yDist + 1e-2;
	double xMax = bound[2] + xDist;
	double yMax = bound[3] + yDist;

	double [] outerBoundary = new double[ 8 ];

	outerBoundary[0] = xMax;	outerBoundary[1] = yMax;
	outerBoundary[2] = xMin;	outerBoundary[3] = yMax;
	outerBoundary[4] = xMin;	outerBoundary[5] = yMin;
	outerBoundary[6] = xMax;	outerBoundary[7] = yMin;

	return outerBoundary;
    }

    public static double [][] getCirclePair( Schottky schottky, int index, int circleDiscr ) {

	final double [][] circlePair = new double[2][ 2*circleDiscr ];

	final Complex z = new Complex();
	final Complex Z = new Complex();
	
	final double r = schottky.getRadius( index );
	     
	final double cx = schottky.center[index][0].re;
	final double cy = schottky.center[index][0].im;

	for( int j=0, k=0; j<circleDiscr; j++ ) {
	    circlePair[0][k++] = cx + r * Math.sin( 2 * j * Math.PI / circleDiscr );
	    circlePair[0][k++] = cy + r * Math.cos( 2 * j * Math.PI / circleDiscr );
	}
	     	
	for( int j=0, k=0; j<circleDiscr; j++ ) {
	    z.assign(circlePair[0][k],
		     circlePair[0][k+1] );  
		
	    schottky.generator[index].applyTo( z, Z );

	    circlePair[1][k++] = Z.re;
	    circlePair[1][k++] = Z.im;		
	}

	return circlePair;
    }

    protected double [] getOuterBoundary() {
    	return getOuterBoundary( schottky, xDist, yDist );
    }
    
    public double [][] getBoundary() {

    	double [][] boundary =null;
    	
	if( boundary == null || schottky.getNumGenerators() * 2 + 1 != boundary.length )
	    boundary = new double[schottky.getNumGenerators() * 2 + 1][];

	boundary[0] = getOuterBoundary();

	for( int i=0, k=1; i < schottky.getNumGenerators(); i++ ) {
	    
	    double [][] circlePair = getCirclePair( schottky, i, circleDiscr );
	    
	    boundary[k++] = circlePair[0];
	    boundary[k++] = circlePair[1];
	}
	return boundary;
    }

    public void updateTriangulator( final boolean rebuildTriangulator ) {
	
	if( triangulator == null || rebuildTriangulator )
	    triangulator = new Ruppert( boundary );

	triangulator.setAngleConstraint( minAngle );
	triangulator.setAreaConstraint( maxArea );
	
	triangulator.refine();

	point = triangulator.getPoints();
    }

    public int getLengthOfBoundaryComponent(int i ) {
    	return boundary[i].length;
    }
    
    public void update() {

	boundary = getBoundary();

	updateTriangulator( true );
    }

    public void addPoint( final double x, final double y ) {
	triangulator.addPoint( x, y );

	System.out.println( "adding ("+x+","+y+")" );

	updateTriangulator( false );
    }

    public void addVerticalLine( final double x, final int resolution ) {

	final double step = ( yMax - yMin ) / resolution;

	double y=yMax;

	for( int i=0; i<resolution; i++, y -= step )
	    addPoint( x, y );

	triangulator.addPoint( x, yMin );

	updateTriangulator( false );
    }

    public void refine() {
	updateTriangulator( false );
    }

    public int getNumOfPoints() {
	return triangulator.getNumPoints();
    }

    public int getNumOfElements() {
	return triangulator.getNumFaces();
    }

    public double [] getPoints() {
	return point;
    }

    public int [] getIndices() {
	return triangulator.getIndices();
    }

    public int [] getNeighbors() {
	return triangulator.getNeighbors();
    }
}      
	




