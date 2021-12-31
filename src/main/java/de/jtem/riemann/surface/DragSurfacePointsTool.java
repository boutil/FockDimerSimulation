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

import java.awt.event.MouseEvent;
import java.util.Vector;

import de.jtem.moebiusViewer.AbstractPicker;
import de.jtem.moebiusViewer.shape.Point;
import de.jtem.moebiusViewer.tool.AbstractTool;
import de.jtem.riemann.constrainedComplex.IllegalCoordinateException;

public class DragSurfacePointsTool extends AbstractTool{
 
    private static final long serialVersionUID = 1L;

    SurfacePoint pickedSurfacePoint;

    Point pickedPoint;

    RamifiedCoveringShape selection;
    
    /**
       * Get the value of pickedSurfacePoint.
       * @return Value of pickedSurfacePoint.
       */
    public SurfacePoint getPickedSurfacePoint() {return pickedSurfacePoint;}
    
    /**
       * Set the value of pickedSurfacePoint.
       * @param v  Value to assign to pickedSurfacePoint.
       */
    public void setPickedSurfacePoint(SurfacePoint  v) {this.pickedSurfacePoint = v;}

    /**
       * Get the value of pickedSurfacePoint.
       * @return Value of pickedSurfacePoint.
       */
    public Point getPickedPoint() {return pickedPoint;}
    
    /**
       * Set the value of pickedPoint.
       * @param v  Value to assign to pickedPoint.
       */
    public void setPickedPoint(Point  v) {this.pickedPoint = v;}
    
    
    public DragSurfacePointsTool() {
	setLabel(  "drag surface points" );
    }
    
    public void mousePressed(MouseEvent e) {
	
	Vector dummy = new Vector();

	super.mousePressed( e );
	
	selection = (RamifiedCoveringShape)(getViewer().getSelection());
	
	int j = 1;
	
	AbstractPicker pickerContext = viewer.getPickerContext();

	pickerContext.setTransformStack( transformPath );

	pickerContext.initPicking ((int) e.getX (), (int) e.getY ());
	
	for( int i=0; i<selection.surface.numOfBranchPoints; i++, j++ ) {
	    
	    pickerContext.point( selection.surface.branchPoint[i].getRe(),
				 selection.surface.branchPoint[i].getIm() );
	    
	    if( pickerContext.isPicked() ) {
		
		pickedSurfacePoint = selection.surface.branchPoint[i];
		
		pickedPoint = selection.point[j]; 

		return;
	    }
	}
	
	for( int i=0; i<selection.surface.numOfSingularPoints; i++, j++ ) {
	    
	    pickerContext.point( selection.surface.singularPoint[i].getRe(),
				 selection.surface.singularPoint[i].getIm() );
	    
	    if( pickerContext.isPicked() ) {
		
		pickedSurfacePoint = selection.surface.singularPoint[i];
		
		pickedPoint = selection.point[j]; 

		return;
	    }
	}

	for( int i=0; i<selection.surface.numOfDistinguishedPoints; i++, j++ ) {
	    
	    pickerContext.point( selection.surface.distinguishedPoint[i].getRe(),
				 selection.surface.distinguishedPoint[i].getIm() );
	    
	    if( pickerContext.isPicked() ) {
		
		pickedSurfacePoint = selection.surface.distinguishedPoint[i];
		
		pickedPoint = selection.point[j]; 

		return;
	    }
	}
	
	pickedSurfacePoint = null;
	pickedPoint        = null;
    }

    public void mouseDragged(MouseEvent e) {

	if( pickedSurfacePoint != null ) {
	    super.mouseDragged( e );

	    try {
	   
		pickedSurfacePoint.setProjectedCoords( newPick );
		
		selection.compute();
	    
		getViewer().repaint();
	    }
	    catch( IllegalCoordinateException exc ) {
		java.awt.Toolkit.getDefaultToolkit().beep();
	    }
	}
    }

    public void mouseReleased(MouseEvent e) {
	selection = null;
    }

}




























