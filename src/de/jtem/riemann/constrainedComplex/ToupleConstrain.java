// JTEM - Java Tools for Experimental Mathematics
// Copyright (C) 2001 JEM-Group
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

package de.jtem.riemann.constrainedComplex;

import de.jtem.mfc.field.Complex;

public abstract class ToupleConstrain extends Constrain {

	private static final long serialVersionUID = 1L;

	ConstrainedComplex A, B;

	final Complex newCoordsOfA = new Complex();
	final Complex newCoordsOfB = new Complex();
	final Complex newCoords = new Complex();

	public ConstrainedComplex getA() {
		return A;
	}

	public ConstrainedComplex getB() {
		return B;
	}

	public ToupleConstrain(ConstrainedComplex aA, ConstrainedComplex aB) {

		A = aA;
		B = aB;

		A.setConstrain(this);
		B.setConstrain(this);
	}

	public boolean effects(ConstrainedComplex aConstrainedComplex) {
		return aConstrainedComplex == A || aConstrainedComplex == B;
	}

	public void deconstrain() {

		ConstrainedComplex tmpA = A;
		ConstrainedComplex tmpB = B;

		A = null;
		B = null;

		tmpA.setConstrain(null);
		tmpB.setConstrain(null);
	}

	public abstract void move(ConstrainedComplex aConstrainedComplex,
			Complex newCoords);

	public int getNumOfParameters(ConstrainedComplex aConstrainedComplex) {
		return aConstrainedComplex == A ? 2 : 0;
	}

	public void setByParameter(double[] p, int offset,
			ConstrainedComplex aConstrainedComplex) {
		if (aConstrainedComplex == A) {
			newCoords.assign(p[offset], p[offset + 1]);
			move(A, newCoords);
		}
	}

	public void getValue(double[] p, int offset,
			ConstrainedComplex aConstrainedComplex) {
		if (aConstrainedComplex == A) {
			p[offset] = A.coords.re;
			p[offset + 1] = A.coords.im;
		}
	}

}
