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

/**
 * @author schmies
 *
 * To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
public class ConstrainedComplex extends de.jtem.mfc.field.ComplexValue implements java.io.Serializable, Cloneable {

		private static final long serialVersionUID = 1L;
		protected Complex coords;
		private final Complex newCoords = new Complex();
		protected Constrain constrain;

		public double getRe() {
			return coords.re;
		}

		public double getIm() {
			return coords.im;
		}

		public final Constrain getConstrain() {
			return constrain;
		}

		public final void deconstrain() {
			if (constrain != null)
				constrain.deconstrain();
		}

		public final void setConstrain(Constrain aConstrain) {
		
			if (aConstrain == null) {
		
				if (constrain == null)
					return;
		
				if (constrain.effects(this))
					throw new IllegalArgumentException("constrain still effects surface point");
		
				constrain = null;
		
			} else if (aConstrain.effects(this)) {
		
				if (constrain != null)
					constrain.deconstrain();
		
				constrain = aConstrain;
		
			} else
				throw new IllegalArgumentException("constrain does not effect surface point");
		}

		public int getNumOfParameters() {
			if (constrain != null)
				return constrain.getNumOfParameters(this);
			else
				return 2;
		}

		public void setByParameter(double[] p, int offset) {
			if (constrain != null)
				constrain.setByParameter(p, offset, this);
			else {
				newCoords.assign(p[offset], p[offset + 1]);
				move(newCoords);
			}
		
		}

		public void getValue(double[] p, int offset) {
			if (constrain != null)
				constrain.getValue(p, offset, this);
			else {
				p[offset] = coords.re;
				p[offset + 1] = coords.re;
			}
		}

		public final Complex getCoords(Complex theCoords) {
			if (theCoords == null)
				return new Complex(coords);
		
			theCoords.assign(coords);
		
			return theCoords;
		}

		public final void setCoords(Complex newCoords) throws IllegalCoordinateException {
		
			if (constrain != null)
				constrain.move(this, newCoords);
			else
				move(newCoords);
		}

		public void setProjectedCoords(Complex newCoords) {
			if (constrain != null)
				constrain.projectedMove(this, newCoords);
			else
				move(newCoords);
		
		}

		protected void move(Complex newCoords) {
		}

		protected void setNewCoords(Complex newCoords) {	
			coords.assign(newCoords);
		}
		
}
