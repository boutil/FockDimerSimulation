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

package de.jtem.riemann.util;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;

import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexConstant;
import de.jtem.mfc.field.ComplexValue;

/**
 * @author schmies
 * 
 *  
 */

public class RadialPathGenerator implements Cloneable, Serializable {

	private static final long serialVersionUID = 1L;

	public static final double EPS = 1e-12;

	private TroublePoint[] troublePoints;

	private TroublePoint[] pointsInTheWay;
	int numOfPointsInTheWay;

	private Complex origin;

	public RadialPathGenerator() {
		this( new Complex() );
	}
	
	public RadialPathGenerator(ComplexValue origin) {
		this.origin = new Complex(origin);
	}
	
	public RadialPathGenerator(
		ComplexValue[] troublePoints,
		ComplexValue origin) {
		this(origin);
		setTroublePoints(troublePoints);
	}

	public void setTroublePoints(ComplexValue[] troublePoints) {
		
		if (this.troublePoints == null
				|| this.troublePoints.length != troublePoints.length) {
			this.troublePoints = new TroublePoint[troublePoints.length];
			
			for (int i = 0; i < troublePoints.length; i++) {
				this.troublePoints[i] = new TroublePoint(troublePoints[i]);
			}
			
		} else {
			for (int i = 0; i < troublePoints.length; i++) {
				this.troublePoints[i].assign(troublePoints[i]);
			}
		}
		
		if( pointsInTheWay == null || pointsInTheWay.length != troublePoints.length)
			pointsInTheWay = new TroublePoint[ troublePoints.length ];
		
		update();
	}


	public Complex getOrigin() {
		return origin;
	}

	private Complex tmp = new Complex();

	public double getMinimalDistanceToTroublePoints(ComplexValue aPoint) {

		tmp.assign(aPoint);
		double min_dist = Double.MAX_VALUE; //tmp.distSqr(origin);

		for (int i = 0; i < troublePoints.length; i++) {

			double dist = tmp.distSqr(troublePoints[i]);

			if (dist < min_dist)
				min_dist = dist;
		}

		return Math.sqrt(min_dist);
	}

	public final void updateDistanceRadii() {

		for (int i = 0; i < troublePoints.length; i++) {

			double radius = getMinimalDistance(troublePoints[i]) / 2;

			troublePoints[i].setRadius( Math.min(radius,
					0.95 * origin.dist( troublePoints[i] ) ) );
		}
	}

	double getMinimalDistance(Complex aPoint) {

		double min_dist =  Double.MAX_VALUE; //aPoint.distSqr(origin);

		for (int i = 0; i < troublePoints.length; i++)
			if (troublePoints[i] != aPoint) {

				double dist = aPoint.distSqr(troublePoints[i]);

				if (dist < min_dist)
					min_dist = dist;
			}

		return Math.sqrt(min_dist);
	}

	void computePath(
		TroublePoint aPoint,
		ComplexVector aPath,
		boolean pathIsLoop) {

		
		LineSeg lineToAPoint = new LineSeg(origin, aPoint);

		LineSeg aLineSeg = new LineSeg();

		numOfPointsInTheWay = 0;

		for (int i = 0; i < troublePoints.length; i++) {

			if (troublePoints[i].distSqr(aPoint) > 1e-20) {

				aLineSeg.set(troublePoints[i].box[1], troublePoints[i].box[3]);

				if (lineToAPoint.intersects(aLineSeg)) {

					pointsInTheWay[numOfPointsInTheWay++] = troublePoints[i];
				}
			}
		}

		aPath.newSize(
			pathIsLoop
				? 7 + 6 * numOfPointsInTheWay
				: 2 + 3 * numOfPointsInTheWay + (aPoint.radius == 0 ? 0 : 1));

		for (int i = 0; i < numOfPointsInTheWay; i++) {

			Complex I = new Complex();

			if (aPoint.isToTheRightOf(pointsInTheWay[i])) {

				aLineSeg.set(
					pointsInTheWay[i].box[0],
					pointsInTheWay[i].box[1]);

				if (lineToAPoint.intersects(aLineSeg, I)) {

					aPath.set(1 + 3 * i, I);
					aPath.set(2 + 3 * i, pointsInTheWay[i].box[1]);

					aLineSeg.set(
						pointsInTheWay[i].box[1],
						pointsInTheWay[i].box[2]);

					lineToAPoint.intersects(aLineSeg, I);

					aPath.set(3 + 3 * i, I);
				}

			} else {

				aLineSeg.set(
					pointsInTheWay[i].box[0],
					pointsInTheWay[i].box[3]);

				if (lineToAPoint.intersects(aLineSeg, I)) {

					aPath.set(1 + 3 * i, I);
					aPath.set(2 + 3 * i, pointsInTheWay[i].box[3]);

					aLineSeg.set(
						pointsInTheWay[i].box[3],
						pointsInTheWay[i].box[2]);

					lineToAPoint.intersects(aLineSeg, I);

					aPath.set(3 + 3 * i, I);
				}
			}
		}

		aPath.set(0, origin);

		if (pathIsLoop) {
			int i;

			// copy the way back
			for (i = 0; i < 1 + 3 * numOfPointsInTheWay; i++)
				aPath.set(6 + 6 * numOfPointsInTheWay - i, aPath.get(i));

			// make loop arount aPoint
			for (int j = 0; j < 4; j++, i++)
				aPath.set(i, aPoint.box[j]);

			aPath.set(i, aPoint.box[0]);
		} else {
			// just assign the last entry to aPoint itself
			int i = 1 + 3 * numOfPointsInTheWay;

			// if the radius is zero we do not insert a point on the last
			// segment
			if (aPoint.radius > 0)
				aPath.set(i++, aPoint.box[0]);

			aPath.set(i, aPoint);
		}
	}

	public ComplexVector getLoopAroundPoint(ComplexValue aPoint) {

		ComplexVector aLoop = new ComplexVector(0);

		computeLoopAroundTroublePoint(aPoint, aLoop);

		return aLoop;
	}

	TroublePoint getTroublePointFor( ComplexValue value ) {
		for( int i=0; i<troublePoints.length; i++)
			if( value == troublePoints[i].delegate )
				return troublePoints[i];
		throw new RuntimeException( "given value is not a trouble point");
	}
	
	public void computeLoopAroundTroublePoint(
		ComplexValue aPoint,
		ComplexVector aLoop) {

		computePath( getTroublePointFor(aPoint), aLoop, true);
	}

	public ComplexVector[] getLoopsAroundTroublePoints( ComplexValue [] point) {

		ComplexVector[] loops = new ComplexVector[point.length];

		computeLoopsAroundTroublePoints(point, loops);

		return loops;
	}

	public void computeLoopsAroundTroublePoints(
		ComplexValue[] point,
		ComplexVector[] loop) {

		int n = point.length;

		if (loop.length != n)
			throw new IllegalArgumentException("array for loops has not the size of array of surface points.");

		for (int i = 0; i < n; i++) {

			if (loop[i] == null)
				loop[i] = new ComplexVector(0);


			if( point[i] == troublePoints[i].delegate )
				computePath(troublePoints[i], loop[i], true);
			else 
				computePath( getTroublePointFor( point[i]), loop[i], true );
			
		}
	}

	public ComplexVector getPathToPoint(ComplexValue aPoint) {

		ComplexVector aPath = new ComplexVector(0);

		computePathToPoint(aPoint, aPath);

		return aPath;
	}

	public void computePathToPoint(
		ComplexValue aPoint,
		ComplexVector aPath) {
		computePath( new TroublePoint( aPoint ), aPath, false);
	}

	public ComplexVector[] getPathsToPoints( ComplexValue[] point) {

		ComplexVector[] loops = new ComplexVector[point.length];

		computePathsToPoints(point, loops);

		return loops;
	}

	public void computePathsToPoints( ComplexValue [] point, ComplexVector[] loop) {

		int n = point.length;

		if (loop.length != n)
			throw new IllegalArgumentException("array for loops has not the size of array of surface points.");

		for (int i = 0; i < n; i++) {

			if (loop[i] == null)
				loop[i] = new ComplexVector(0);

			computePathToPoint(point[i], loop[i]);
		}
	}

	public double getDistRadiusOfTroublePoint( ComplexValue value ) {
		return this.getTroublePointFor(value).radius;
	}
	
	public void computeBoxArroundTroublePoint( ComplexValue value, ComplexVector box ) {
		box.newSize(5);
		TroublePoint p = getTroublePointFor(value);
		for( int i=0; i<4; i++ )
			box.set(i,p.box[i]);
		box.set(4,p.box[0]);
	}
	
	public ComplexVector getBoxArroundTroublePoint( ComplexValue value ) {
		ComplexVector box= new ComplexVector(5);
		computeBoxArroundTroublePoint( value, box );
		return box;
	}
	
	class RadialComparator implements Comparator {

		Complex origin, tmp = new Complex();

		RadialComparator(Complex anOrigin) {
			origin = anOrigin;
		}

		public int compare(Object o1, Object o2) {

			tmp.assignMinus(((TroublePoint) o1), origin);
			double a = tmp.absSqr();

			tmp.assignMinus(((TroublePoint) o2), origin);
			double b = tmp.absSqr();

			if (a < b)
				return -1;
			else if (a > b)
				return 1;

			return 0;
		}
	}

	public void update() {

		for (int i = 0; i < troublePoints.length; i++) {
			this.troublePoints[i].assign( troublePoints[i].delegate );
		}
		
		updateDistanceRadii();

		Arrays.sort(troublePoints, new RadialComparator(origin));
	}

	class TroublePoint extends Complex implements Serializable {

		private static final long serialVersionUID = 1L;

		double arg, abs;

		final private Complex tmp = new Complex();
		final private Complex newCoords = new Complex();

		private ComplexValue delegate;

		double radius;

		Complex box[] = new Complex[4];

		TroublePoint(ComplexValue z) {
			super(z);
			delegate = z;

			for (int i = 0; i < 4; i++)
				box[i] = new Complex();
			
			arg = getArgFor(this);
			abs = getAbsFor(this);
			
			computeBox();
		}

		void assign(ComplexValue z) {
			super.assign(z);

			delegate = z;

			arg = getArgFor(this);
			abs = getAbsFor(this);
		}
		
		public final double getDeltaArgTo(TroublePoint otherPoint) {

			double newArg = otherPoint.arg;

			if (newArg < arg)
				newArg += 2 * Math.PI;

			double deltaArg = newArg - arg;

			if (deltaArg > Math.PI)
				deltaArg -= 2 * Math.PI;

			return deltaArg;
		}

		public boolean isToTheRightOf(TroublePoint other) {

			double deltaArg = getDeltaArgTo(other);

			if (deltaArg == 0)
				return abs < other.abs;

			return deltaArg > 0;
		}

		final double getArgFor(Complex someCoords) {
			tmp.assignMinus(someCoords, origin);
			double arg = ComplexValue.argPositive(tmp);
			if (arg > Math.PI)
				arg -= 2 * Math.PI;
			return arg;
		}

		final double getAbsFor(Complex someCoords) {
			tmp.assignMinus(someCoords, origin);
			return tmp.abs();
		}

		private Complex radiator = new Complex();

		void setRadius(double aRadius) {
			radius = aRadius;

			computeBox();
		}

		void computeBox() {

			radiator.assignMinus(origin, this);
			radiator.assignTimes(radius / radiator.abs());

			for (int j = 0;
				j < 4;
				j++, radiator.assignTimes(ComplexConstant.I))
				box[j].assignPlus(this, radiator);
		}

	}

}
