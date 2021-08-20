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

package de.jtem.riemann.surface.hyperElliptic;

import java.awt.Graphics;
import java.beans.Customizer;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;

import de.jtem.moebiusViewer.moebius.MoebiusViewer;
import de.jtem.riemann.surface.RamifiedCoveringShape;

public class HyperEllipticCurveCustomizer extends MoebiusViewer 
  implements Customizer {
  private static final long     serialVersionUID = 1L;
  PropertyChangeSupport         propChange = new PropertyChangeSupport (this);

  public HyperEllipticCurveCustomizer (){}

    RamifiedCoveringShape rps = new RamifiedCoveringShape();

  Object        theObject;
  boolean       firstPaint;

  public void setObject (Object o) {
    if (o != theObject) {
      theObject = o;

      rps.set( (HyperEllipticCurve)o );

      removeAll();
      add (rps);

      firstPaint = true;
    } 
  } 

  public void paint (Graphics g) {
    if (firstPaint) {
      encompass ();

      firstPaint = false;
    } 

    super.paint (g);
  } 

  public void addPropertyChangeListener (PropertyChangeListener listener) {
    propChange.addPropertyChangeListener (listener);
  } 

  public void removePropertyChangeListener (PropertyChangeListener listener) {
    propChange.removePropertyChangeListener (listener);
  } 

  protected void firePropertyChange (java.lang.String s, java.lang.Object o, 
          java.lang.Object n) {
    propChange.firePropertyChange (s, o, n);
    System.out.println ("fire");
  } 

}

