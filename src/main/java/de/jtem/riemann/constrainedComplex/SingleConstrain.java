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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

package de.jtem.riemann.constrainedComplex;

import de.jtem.mfc.field.Complex;

public abstract class SingleConstrain extends Constrain {

	private static final long serialVersionUID = 1L;

	ConstrainedComplex A;

	public ConstrainedComplex getA() {
		return A;
	}

	public boolean effects(ConstrainedComplex aConstrainedComplex) {
		return aConstrainedComplex == A;
	}

	public void deconstrain() {

		ConstrainedComplex tmpA = A;

		A = null;

		tmpA.setConstrain(null);
	}

	public void projectedMove(
		ConstrainedComplex aConstrainedComplex,
		Complex newCoords) {
	}

	public void move(
		ConstrainedComplex aConstrainedComplex,
		Complex newCoords) {

		if (!A.coords.equals(newCoords))
			throw new IllegalCoordinateException(aConstrainedComplex);
	}
}
