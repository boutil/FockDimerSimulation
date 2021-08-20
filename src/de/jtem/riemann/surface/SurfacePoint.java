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

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.field.ComplexValue;
import de.jtem.riemann.constrainedComplex.ConstrainedComplex;

public abstract class SurfacePoint
	extends ConstrainedComplex
	implements Serializable, Cloneable, Comparable {

	private static final long serialVersionUID = 1L;

	int index;

	double arg, abs;

	final private Complex tmp = new Complex();
	private RamifiedCovering surface;

	private Complex origin;

	SurfacePoint(double x, double y) {

		coords = new Complex(x, y);

	}

	SurfacePoint(Complex z) {

		this(z.re, z.im);
	}

	SurfacePoint(RamifiedCovering aSurface, double x, double y) {

		this(x, y);

		setSurface(aSurface);
	}

	public final int getIndex() {
		return index;
	}

	public final Complex getCoords() {
		return new Complex(coords);
	}

	protected void move(Complex newCoords) {
		super.move( newCoords );
		surface.move(this, newCoords);
	}

	protected void setNewCoords(Complex newCoords) {
		super.setNewCoords( newCoords );

		arg = getArgFor(coords);
		abs = getAbsFor(coords);
	}

	public final RamifiedCovering getSurface() {
		return surface;
	}

	public final double getArgFor(Complex someCoords) {
		tmp.assignMinus(someCoords, origin);
		double arg = ComplexValue.argPositive(tmp);
		if (arg > Math.PI)
			arg -= 2 * Math.PI;
		return arg;
	}

	public final double getAbsFor(Complex someCoords) {
		tmp.assignMinus(someCoords, origin);
		return tmp.abs();
	}

	public final double getArg() {
		return arg;
	}

	public final double getAbs() {
		return abs;
	}

	public void setSurface(RamifiedCovering aSurface) {

		if (surface != null)
			throw new SecurityException("the surface can only be set once");

		surface = aSurface;

		origin = surface.getOrigin();

		arg = getArgFor(coords);
		abs = getAbsFor(coords);
	}

	public final double getDeltaArgTo(Complex newCoords) {

		double newArg = getArgFor(newCoords);

		if (newArg < arg)
			newArg += 2 * Math.PI;

		double deltaArg = newArg - arg;

		if (deltaArg > Math.PI)
			deltaArg -= 2 * Math.PI;

		return deltaArg;
	}

	public final double getDeltaArgTo(SurfacePoint other) {

		double newArg = other.arg;

		if (newArg < arg)
			newArg += 2 * Math.PI;

		double deltaArg = newArg - arg;

		if (deltaArg > Math.PI)
			deltaArg -= 2 * Math.PI;
		return deltaArg;
	}

	public boolean isToTheRightOf(Complex newCoords) {

		double deltaArg = getDeltaArgTo(newCoords);

		if (deltaArg == 0)
			return abs < getAbsFor(newCoords);

		return deltaArg > 0;
	}

	public boolean isToTheRightOf(SurfacePoint other) {

		double deltaArg = getDeltaArgTo(other);

		if (deltaArg == 0)
			return abs < other.abs;

		return deltaArg > 0;
	}

	public int compareTo(Object aPoint) {

		SurfacePoint other = (SurfacePoint) aPoint;

		if (arg < other.arg)
			return -1;
		else if (arg > other.arg)
			return 1;
		else if (abs < other.abs)
			return -1;
		else if (abs > other.abs)
			return 1;

		return 0;
	}

	public static final int getIndexOfSurfacePoint(
		SurfacePoint[] somePoints,
		SurfacePoint aPoint) {

		for (int i = 0; i < somePoints.length; i++)
			if (somePoints[i] == aPoint)
				return i;

		return -1;
	}

	public static final int getIndexOfSurfacePointWithCoords(
		SurfacePoint[] somePoints,
		Complex coords,
		double epsSqr) {

		for (int i = 0; i < somePoints.length; i++) {
			if (coords.distSqr(somePoints[i].coords) < epsSqr)
				return i;
		}

		return -1;
	}

	public static final int getIndexOfSurfacePointWithCoords(
		SurfacePoint[] somePoints,
		Complex coords) {

		return getIndexOfSurfacePointWithCoords(
			somePoints,
			coords,
			RamifiedCovering.EPS * RamifiedCovering.EPS);
	}

	public static final SurfacePoint getSurfacePointWithCoords(
		SurfacePoint[] somePoints,
		Complex coords,
		double epsSqr) {

		for (int i = 0; i < somePoints.length; i++) {
			if (coords.distSqr(somePoints[i].coords) < epsSqr)
				return somePoints[i];
		}

		return null;
	}

	public static SurfacePoint[] getSurfacePoints(Complex[] specialPoints) {
		SurfacePoint[] somePoints = new SurfacePoint[specialPoints.length];

		for (int i = 0; i < somePoints.length; i++)
			somePoints[i] = new SingularPoint(specialPoints[i]);

		return somePoints;
	}

}
