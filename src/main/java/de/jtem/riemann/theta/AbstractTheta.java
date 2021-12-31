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

package de.jtem.riemann.theta;

import java.io.Serializable;

import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;

/** Abstract super class for all theta function derivatives. */
abstract class AbstractTheta implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    private  final Complex tmp             = new Complex();

    private  final Complex factor          = new Complex();

    private  final Complex exponent        = new Complex();
    private  final Complex expOfExponent   = new Complex();

    private  final Complex thetaSumZ       = new Complex();
    private  final Complex thetaSumX       = new Complex();
    private  final Complex thetaSumY       = new Complex();
    private  final Complex thetaSumXY      = new Complex();

    private  final Complex i2Pi            = new Complex( 0, 2*Math.PI);

    abstract public void theta( final ComplexVector Z, final Complex factor, final Complex thetaSumZ );

    /**
     * Evaluates the Riemann theta function at <code>Z</code>.
     * <p align=center>
     *   <code>
     *     &theta;(z|B) = thetaAtZ
     *   </code>
     * </p>
     * @param Z argument vector
     * @param thetaAtZ value of the Riemann theta function at <code>Z</code>
     */
    public final void theta( final ComplexVector Z, final Complex thetaAtZ ) {

	theta( Z, exponent, thetaSumZ );

        expOfExponent.assignExp( exponent );
	thetaAtZ.assignTimes( thetaSumZ, expOfExponent );
    }


    /**
     * Returns the Riemann theta function at <code>Z</code>.
     * @param Z argument vector
     * @return value of the Riemann theta function at <code>Z</code>
     */
    public final Complex theta( final ComplexVector Z ) {

	theta( Z, exponent, thetaSumZ );

        expOfExponent.assignExp( exponent );
	return thetaSumZ.times( expOfExponent );
    }

    abstract public void dTheta( final ComplexVector Z,
				 final ComplexVector Y,
				 final Complex factor,
				 final Complex thetaSumZ,
				 final Complex thetaSumY );

    /**
     * Evaluates the Riemann theta function and its first derivative
     * in the <code>X</code> direction at <code>Z</code>.
     * <p align=center>
     * <code>
     * <table>
     * <tr> <td>              &theta;(z|B) </td> <td> = </td> <td> thetaAtZ  </td> </tr>
     * <tr> <td> D<sub>X</sub>&theta;(z|B) </td> <td> = </td> <td> DXThetaAtZ</td> </tr>
     * </table>
     * </code>
     * <p>
     * Evaluating the Riemann theta funciton and its first derivative simultanesly
     * can be performed with almost no extra cost.
     * @param Z argument vector
     * @param X direction of derivative
     * @param thetaAtZ value of the Riemann theta function at <code>Z</code>
     * @param DXThetaAtZ value of first derivative in <code>X</code> direction
     */
    public final void dTheta( final ComplexVector Z,
			      final ComplexVector X,
			      final Complex thetaAtZ,
			      final Complex DXThetaAtZ ) {

	dTheta( Z, X, exponent, thetaAtZ, DXThetaAtZ );

        expOfExponent.assignExp( exponent );
	thetaAtZ.assignTimes( expOfExponent );
	DXThetaAtZ.assignTimes( expOfExponent );
    }

    /**
     * Returns the first derivative of the Riemann theta function
     * in the <code>X</code> direction at <code>Z</code>.
     * @param Z argument vector
     * @param X direction of derivative
     * @return value of first derivative in <code>X</code> direction
     */
    public final Complex dTheta( final ComplexVector Z, final ComplexVector X ) {
	dTheta( Z, X, exponent, thetaSumZ, thetaSumX );

        expOfExponent.assignExp( exponent );
        return thetaSumX.times( expOfExponent );
    }

    abstract public void ddTheta( final ComplexVector Z,
				  final ComplexVector Y,
				  final ComplexVector X,
				  final Complex factor,
				  final Complex thetaSumZ,
				  final Complex thetaSumX,
				  final Complex thetaSumY,
				  final Complex thetaSumXY );

    /**
     * Evaluates the Riemann theta function, its first derivative
     * in the <code>X</code> and <code>Y</code> direction and,
     * its second derivative into the same direction at <code>Z</code>.
     * <p align=center>
     * <code>
     * <table>
     * <tr> <td>                 &theta;(z|B)      </td> <td> = </td> <td> thetaAtZ    </td> </tr>
     * <tr> <td> D<sub>X</sub>   &theta;(z|B)      </td> <td> = </td> <td> DXThetaAtZ  </td> </tr>
     * <tr> <td> D<sub>Y</sub>   &theta;(z|B)      </td> <td> = </td> <td> DYThetaAtZ  </td> </tr>
     * <tr> <td> D&sup2;<sub>X,Y</sub>&theta;(z|B) </td> <td> = </td> <td> DXDYThetaAtZ</td> </tr>
     * </table>
     * </code>
     * <p>
     * Evaluating the Riemann theta funciton and its first derivative simultanesly
     * can be performed with almost no extra cost.
     * @param Z argument vector
     * @param X direction of derivative
     * @param Y direction of derivative
     * @param thetaAtZ oscillatory part
     * @param DXThetaAtZ oscillatory part of first derivative in <code>X</code> direction
     * @param DYThetaAtZ oscillatory part of first derivative in <code>Y</code> direction
     * @param DXDYThetaAtZ oscillatory part of second derivative
     */
    public final void ddTheta( final ComplexVector Z,
			       final ComplexVector X,
			       final ComplexVector Y,
			       final Complex thetaAtZ,
			       final Complex DXThetaAtZ,
			       final Complex DYThetaAtZ,
			       final Complex DXDYThetaAtZ ) {

	ddTheta( Z, X, Y, exponent, thetaAtZ, DXThetaAtZ, DYThetaAtZ, DXDYThetaAtZ );

        expOfExponent.assignExp( exponent );
	thetaAtZ.assignTimes( expOfExponent );
	DXThetaAtZ.assignTimes(  expOfExponent );
	DYThetaAtZ.assignTimes(  expOfExponent );
	DXDYThetaAtZ.assignTimes( expOfExponent );
    }

    /**
     * Evaluates the Riemann theta function, its first derivative
     * in the <code>X</code> and <code>Y</code> direction and,
     * its second derivative into the same direction at <code>Z</code>.
     * @param Z argument vector
     * @param X direction of derivative
     * @param Y direction of derivative
     * @param thetaAtZ oscillatory part
     * @param DXThetaAtZ oscillatory part of first derivative in <code>X</code> direction
     * @param DYThetaAtZ oscillatory part of first derivative in <code>Y</code> direction
     * @param DXDYThetaAtZ oscillatory part of second derivative
     */
    public final Complex ddTheta( final ComplexVector Z, final ComplexVector X, final ComplexVector Y ) {

	ddTheta( Z, X, Y, exponent, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY );

	expOfExponent.assignExp( exponent );

	return thetaSumXY.times( expOfExponent );
    }

    /**
     * Evaluates the first logarythmic derivative of the Riemann theta funciton
     * in the <code>X</code> direction at <code>Z</code>.
     * <p align=center>
     *   <code>
     *     D<sub>X</sub>log( &theta;(z|B) ) = DXLogThetaAtZ
     *   </code>
     * </p>
     * @param Z argument vector
     * @param X direction of derivative
     * @param DXLogThetaAtZ first logarythmic derivative in <code>X</code> direction
     */
    public final void dLogTheta( final ComplexVector Z,
				 final ComplexVector X,
				 final Complex DXLogThetaAtZ ) {

	dTheta( Z, X, exponent, thetaSumZ, thetaSumX );

	DXLogThetaAtZ.assignDivide( thetaSumX, thetaSumZ );
    }

    /**
     * Returns the first logarythmic derivative of the Riemann theta function
     * in the <code>X</code> direction at <code>Z</code>.
     * @param Z argument vector
     * @param X direction of derivative
     * @return first logarythmic derivative in <code>X</code> direction
     */
    public final Complex dLogTheta( final ComplexVector Z, final ComplexVector X ) {
	dTheta( Z, X, factor, thetaSumZ, thetaSumX );

	return thetaSumX.divide( thetaSumZ );
    }

    /**
     * Evaluates the second logarythmic derivative of the Riemann theta funciton
     * in the <code>X</code> and the <code>Y</code> direction at <code>Z</code>.
     * <p align=center>
     *   <code>
     *     D&sup2;<sub>X,Y</sub>log( &theta;(z|B) ) = DXDYLogThetaAtZ
     *   </code>
     * </p>
     * @param Z argument vector
     * @param X direction of derivative
     * @param Y direction of derivative
     * @param DXDYLogThetaAtZ second logarythmic derivative in the <code>X</code> and the <code>Y</code> direction
     */
    public final void ddLogTheta( final ComplexVector Z,
				  final ComplexVector X,
				  final ComplexVector Y,
				  final Complex DXDYLogThetaAtZ ) {

	ddTheta( Z, X, Y, exponent, thetaSumZ, thetaSumX, thetaSumY, thetaSumXY );

	DXDYLogThetaAtZ.assignDivide( thetaSumXY, thetaSumZ );

	tmp.assignTimes( thetaSumX, thetaSumY );
	tmp.assignDivide( thetaSumZ );
	tmp.assignDivide( thetaSumZ );

	DXDYLogThetaAtZ.assignMinus( tmp );
    }

    /**
     * Evaluates the first logarythmic derivative of the Riemann theta funciton
     * in the <code>X</code> and <code>Y</code> direction at <code>Z</code>.
     * <p align=center>
     *   <code>
     *     D&sup2;<sub>X,Y</sub>log( &theta;(z|B) ) = DXDYLogThetaAtZ
     *   </code>
     * </p>
     * @param Z argument vector
     * @param X direction of derivative
     * @param Y direction of derivative
     * @return second logarythmic derivative in <code>X</code> and <code>Y</code> direction
     */
    public final Complex ddLogTheta( final ComplexVector Z,
				     final ComplexVector X,
				     final ComplexVector Y ) {

	final Complex DXDYLogThetaAtZ = new Complex();

	ddLogTheta( Z, X, Y, DXDYLogThetaAtZ );

	return DXDYLogThetaAtZ;
    }
}
















