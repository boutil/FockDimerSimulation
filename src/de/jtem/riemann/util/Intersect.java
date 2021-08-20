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

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.io.Serializable;

public final class Intersect implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;
 
    private static double EPS = 1e-14;

  
  static double lines( Line2D.Double L1, Line2D.Double L2, Point2D.Double I, Point2D.Double coords ) {

    double U0 = L2.x1 - L1.x1;
    double U1 = L2.y1 - L1.y1;

    double V0 = L1.x2 - L1.x1;
    double V1 = L1.y2 - L1.y1;

    double W0 = L2.x1 - L2.x2;
    double W1 = L2.y1 - L2.y2;
  
    double det = V0*W1 - V1*W0; 

    if( det < EPS && det > - EPS )
      return 0;

    if( coords != null  || I != null ) {

      double coords0 = ( W1*U0 - W0*U1 ) / det;
      double coords1 = (-V1*U0 + V0*U1 ) / det;
      
      if( I != null ) {
	I.x = coords0 * L1.x2 + (1 - coords0) * L1.x1;
	I.y = coords0 * L1.y2 + (1 - coords0) * L1.y1;
      }

      if( coords != null ) {
	coords.x = coords0;
	coords.y = coords1;
      }
    }

    return -det;
  }

  static double lineSegments( Line2D.Double L1, Line2D.Double L2, Point2D.Double I, Point2D.Double coords ) {

      if( coords == null )
	  coords = new Point2D.Double();

      double det = lines( L1, L2, I, coords );

      if( det != 0
	  && coords.x < 1+EPS && coords.x > -EPS
	  && coords.y < 1+EPS && coords.y > -EPS ) return det;

      return 0;
  }

}
