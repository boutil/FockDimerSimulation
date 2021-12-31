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

package de.jtem.riemann.surface.homologie;

import java.io.Serializable;

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.blas.IntegerMatrix;
import de.jtem.blas.IntegerVector;
import de.jtem.mfc.field.Complex;
/**
    All members of this class are defined during instantiation</br>
    and are recalculated during the instantiation of
    <i>SymmetricBasis:</i></br>
    </br>
    <i>initialACycles</i></br>
    <i>initialBCycles</i></br>

*/

public class CanonicalBasis extends Basis implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    int genus;

    IntegerMatrix initialACycles;
    IntegerMatrix initialBCycles;


    public void assign( CanonicalBasis aTransform ) {
	super.assign( aTransform );

	genus = aTransform.genus;

	if( aTransform.initialACycles == null ) initialACycles = null;
	else if( initialACycles == null )
	    initialACycles = new IntegerMatrix ( aTransform.initialACycles  );
	else initialACycles.assign(  aTransform.initialACycles  );

	if( aTransform.initialBCycles == null ) initialBCycles = null;
	else if( initialBCycles == null )
	    initialBCycles = new IntegerMatrix ( aTransform.initialBCycles  );
	else initialBCycles.assign(  aTransform.initialBCycles  );

    }

    public Object clone() {
	return new CanonicalBasis( this );
    }

    /**
     * Get the value of genus.
     * @return Value of genus.
     */
    public int getGenus() {return genus;}

    /**
     * Get the value of initialACycles.
     * @return Value of initialACycles.
     */
    public IntegerMatrix getInitialACycles() {return initialACycles;}

    /**
     * Get the value of initialBCycles.
     * @return Value of initialBCycles.
     */
    public IntegerMatrix getInitialBCycles() {return initialBCycles;}

    static IntegerMatrix transformToCanonicForm ( IntegerMatrix scal_prod ) {

	int i, j, k, k1, k2, k3;

	int numOfCycles = scal_prod.getNumRows();

	IntegerMatrix can_basis = IntegerMatrix.id( numOfCycles );

	int zero_space_dim=0;
	for(k=0; k< numOfCycles-1; k+=2) {

	    while( true ) {
		int [] pos = new int[2];
		int [] val = new int[2];

		GuIntegerMatrix.find2LargestInRowStartingFrom(scal_prod, k, k+1, pos, val );

		if(pos[0] == -1) { /* We move the zero line to the end of the matrix */
		    zero_space_dim++;
		    if(zero_space_dim+k>=numOfCycles) break;
		    GuIntegerMatrix.moveColAndRowToEnd (scal_prod, k );
		    GuIntegerMatrix.moveRowToEnd          (can_basis, k );

		} else {

		    if(pos[1] != -1) {

			if(Math.abs(val[0])==Math.abs(val[1])&&pos[0]==k+1)
			    { pos[0]=pos[1]; pos[1]=k+1; k3=val[0]; val[0]=val[1]; val[1]=k3; }
			k3=-val[0]/val[1];
			GuIntegerMatrix.addColAndRow    (scal_prod, pos[0], pos[1], k3 );
			GuIntegerMatrix.addRow             (can_basis, pos[0], pos[1], k3 );

		    } else {

			if(pos[0] != k+1) {

			    GuIntegerMatrix.swapColsAndRows(scal_prod, k+1, pos[0]);
			    GuIntegerMatrix.swapRows          (can_basis, k+1, pos[0]);

			} else {
			    /* We procced the line number k+1 */
			    GuIntegerMatrix.find2LargestInRowStartingFrom (scal_prod, k+1, k+1, pos, val );
			    if( pos[0] == -1 )
				break;

			    while(true) {

				GuIntegerMatrix.find2LargestInRowStartingFrom (scal_prod, k+1, k, pos, val );

				if( pos[1] == -1 && pos[0]==k )
				    break;

				if( pos[1] != -1 ) {

				    if(Math.abs(val[0])==Math.abs(val[1])&&pos[0]==k)
					{pos[0]=pos[1]; pos[1]=k;k3=val[0]; val[0]=val[1]; val[1]=k3;}
				    k3=-val[0]/val[1];
				    GuIntegerMatrix.addColAndRow    (scal_prod, pos[0], pos[1], k3 );
				    GuIntegerMatrix.addRow             (can_basis, pos[0], pos[1], k3 );

				} else {
				    GuIntegerMatrix.swapColsAndRows(scal_prod, k, pos[0]);
				    GuIntegerMatrix.swapRows          (can_basis, k, pos[0]);
				}
			    }
			}
		    }
		}
	    }
	}

	for(k=0; k< numOfCycles-1; k+=2) {
	    if (scal_prod.re[k][k+1] < 0 ) {
		GuIntegerMatrix.multColAndRow(scal_prod, k, -1 );
		GuIntegerMatrix.multRow      (can_basis, k, -1 );
	    }
	}

	return can_basis;
    }

    protected CanonicalBasis() {
    }

    public CanonicalBasis( CanonicalBasis aCanonicalBasis ) {
	this.assign( aCanonicalBasis );
    }

    public CanonicalBasis( IntegerMatrix monodrom, IntegerMatrix singular, IntegerVector disting ) {

	    super( monodrom, singular, disting );

	    /* Now we transform the intersection matrix to the canonic form */

	    if( debug )
		firstCycleIntersections.print("Intersection Matrix");

	    IntegerMatrix can_basis = transformToCanonicForm ( new IntegerMatrix ( firstCycleIntersections ) );

	    if( debug )
		can_basis.print("canonic basis");

	    /* A test -- we calculate the intersection in new basis */
	    /* It can be removed when not needed It is used for debugging purposes only */

	    if( debug ) {
		IntegerMatrix scal_prod
		    = can_basis.times( firstCycleIntersections ).times( can_basis.transpose() );

		scal_prod.print( " Intersection matrix in canonical basis ");
	    }

	    /* end of test */

	    genus = numOfCycles/2;

	    initialACycles = new IntegerMatrix ( genus, numOfEdges );
	    initialBCycles = new IntegerMatrix ( genus, numOfEdges );

	    /*  Now we calculete the canonical basis !!!!!!!!!!!!!!!!!!!!!! */

	    for(int i=0;i<genus;i++)
		for(int j=0;j<numOfEdges;j++)
		    for(int k=0;k<2*genus;k++) {
			initialACycles.re[i][j] += can_basis.re[2*i  ][k] * firstCycleBasis.re[k][j];
			initialBCycles.re[i][j] += can_basis.re[2*i+1][k] * firstCycleBasis.re[k][j];
		    }

	    if( debug ) {
		initialACycles.print( " Initial A-cycles " );
		initialBCycles.print( " Initial B-cycles " );
	    }
    }


    public ComplexMatrix getAPeriodsForEdgeIntegrals( ComplexMatrix I ) {

	ComplexMatrix aPeriod = new ComplexMatrix( genus );

	computeAPeriodsForEdgeIntegrals( I, aPeriod );

	return aPeriod;
    }

    private ComplexMatrix integralMatrix = new ComplexMatrix();
    private ComplexVector integralVector = new ComplexVector();

    public void computeAPeriodsForEdgeIntegrals( ComplexMatrix I, ComplexMatrix aPeriod ) {

      I.getBlock(0,0, numOfEdges, genus, integralMatrix );

      aPeriod.assignTimes( initialACycles, integralMatrix );

      aPeriod.assignTranspose();
    }

    public ComplexVector getAPeriodsForEdgeIntegrals( ComplexVector I ) {

	ComplexVector aPeriod = new ComplexVector( genus );

	computeAPeriodsForEdgeIntegrals( I, aPeriod );

	return aPeriod;
    }

    public void computeAPeriodsForEdgeIntegrals( ComplexVector I, ComplexVector aPeriod ) {

        I.getBlock( 0, numOfEdges, integralVector );

        aPeriod.assignTimes( initialACycles, integralVector );
    }


    public ComplexMatrix getBPeriodsForEdgeIntegrals( ComplexMatrix I ) {

	ComplexMatrix bPeriod = new ComplexMatrix( genus );

	computeBPeriodsForEdgeIntegrals( I, bPeriod );

	return bPeriod;
    }

    public void computeBPeriodsForEdgeIntegrals( ComplexMatrix I, ComplexMatrix bPeriod ) {

      I.getBlock(0, 0, numOfEdges, genus, integralMatrix);

      bPeriod.assignTimes(initialBCycles, integralMatrix);

      bPeriod.assignTranspose();
    }

    public ComplexVector getBPeriodsForEdgeIntegrals( ComplexVector I ) {

	ComplexVector bPeriod = new ComplexVector( genus );

	computeBPeriodsForEdgeIntegrals( I, bPeriod );

	return bPeriod;
    }

    public void computeBPeriodsForEdgeIntegrals( ComplexVector I, ComplexVector bPeriod ) {

      I.getBlock(0, numOfEdges, integralVector);

      bPeriod.assignTimes( initialBCycles, integralVector );

    }


    /*
    public ComplexMatrix getBPeriodsMatrixForEdgeIntegrals( ComplexMatrix I ) {

	ComplexMatrix periodMatrix = new ComplexMatrix( genus, genus );

	computePeriodMatrixForEdgeIntegrals( I, periodMatrix );

	return periodMatrix;
    }
    */

    public void computePeriodMatrixForEdgeIntegrals( ComplexMatrix I, ComplexMatrix periodMatrix ) {

	ComplexMatrix aPeriod = getAPeriodsForEdgeIntegrals( I );
	ComplexMatrix bPeriod = getBPeriodsForEdgeIntegrals( I );

	ComplexMatrix formTransform = aPeriod.invert();

	formTransform.assignTimes( new Complex( 0, 2*Math.PI ) );

	periodMatrix.assignTimes( formTransform, bPeriod );
    }


}










