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

import de.jtem.blas.IntegerMatrix;

final class GuIntegerMatrix implements Serializable, Cloneable {

    static final void swapCols( IntegerMatrix matrix, int I, int J ) {
	matrix.checkColIndex( I );
	matrix.checkColIndex( J );


	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int k=0;k<numOfRows;k++)
	    {
		int buf=matrix.re[k][I];

		matrix.re[k][I]=matrix.re[k][J];
		matrix.re[k][J]=buf;
	    }
    }

    static final void swapRows( IntegerMatrix matrix, int I, int J ) {
	matrix.checkRowIndex( I );
	matrix.checkRowIndex( J );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int k=0;k<numOfCols;k++)
	    {
		int buf=matrix.re[I][k];

		matrix.re[I][k]=matrix.re[J][k];
		matrix.re[J][k]=buf;
	    }
    }

    static final void swapColsAndRows( IntegerMatrix matrix, int I, int J ) {
	swapCols( matrix, I, J );
	swapRows(    matrix, I, J );
    }


    static final void moveColToEnd( IntegerMatrix matrix, int I ) {
	matrix.checkColIndex( I );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int k=0;k<numOfRows;k++)
	    {
		int buf=matrix.re[k][I];

		for(int l=I;l<numOfCols-1;l++)
		    matrix.re[k][l]=matrix.re[k][l+1];

		matrix.re[k][numOfCols-1]=buf;
	    }
    }

    static final void moveRowToEnd( IntegerMatrix matrix, int I ) {
	matrix.checkRowIndex( I );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int k=0;k<numOfCols;k++)
	    {
		int buf=matrix.re[I][k];

		for(int l=I;l<numOfRows-1;l++)
		    matrix.re[l][k]=matrix.re[(l+1)][k];

		matrix.re[(numOfRows-1)][k]=buf;
	    }
    }

    static final void moveColAndRowToEnd ( IntegerMatrix matrix, int I ) {
	moveColToEnd( matrix, I );
	moveRowToEnd( matrix, I );
    }

    static final void addCol( IntegerMatrix matrix, int I, int J, int K ) {
	/* The program adds col J times K to the col I   */
	matrix.checkColIndex( I );
	matrix.checkColIndex( J );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int i1=0;i1<numOfRows;i1++)
	    matrix.re[i1][I]+=(K*matrix.re[i1][J]);
    }

    static final void addRow( IntegerMatrix matrix, int I, int J, int K ) {
	/* The program adds row J times K to the row I   */
	matrix.checkRowIndex( I );
	matrix.checkRowIndex( J );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int i1=0;i1<numOfCols;i1++)
	    matrix.re[I][i1]+=(K*matrix.re[J][i1]);
    }

    static final void addColAndRow( IntegerMatrix matrix, int I, int J, int K ) {
	/* The program adds col J times K to the col I  and the row J times K to the row I */
	addCol( matrix, I, J, K );
	addRow(    matrix, I, J, K );
    }

    static final void multCol          ( IntegerMatrix matrix, int I, int K ) {
	/* The program Multipliess col I to number K  */
	matrix.checkColIndex( I );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int i1=0;i1<numOfRows;i1++)
	    matrix.re[i1][I]*=K;
    }

    static final void multRow             ( IntegerMatrix matrix, int I, int K ) {
	/* The program Multipliess row I to number K  */
	matrix.checkRowIndex( I );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int i1=0;i1<numOfCols;i1++)
	    matrix.re[I][i1]*=K;
    }

    static final void multColAndRow    ( IntegerMatrix matrix, int I, int K ) {
	/* The program Multipliess col I to number  K and the row I to K */
	multCol( matrix, I, K );
	multRow( matrix, I, K );
    }


    static final void divideCol          ( IntegerMatrix matrix, int I, int K ) {
	/* The program divides col I to number K  in a stupid way (witghout checking the rest) */
	matrix.checkColIndex( I );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int i1=0;i1<numOfRows;i1++)
	    matrix.re[i1][I]/=K;
    }



    static final void find2LargestInColStartingFrom (IntegerMatrix matrix, int col, int startPoint, int[] pos, int[] val ) {
	pos[0]=-1;pos[1]=-1;val[0]=0;val[1]=0;

	matrix.checkColIndex( col );
	matrix.checkRowIndex( startPoint );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int j=startPoint; j<numOfRows; j++)
	    {
		int k3=matrix.re[j][col];
		if( k3 != 0 )
		    {
			if( Math.abs(k3) > Math.abs(val[0]) )
			    {
				val[1]=val[0];pos[1]=pos[0];val[0]=k3;pos[0]=j;
			    } else if( Math.abs(k3) > Math.abs(val[1]) ) {  val[1]=k3; pos[1]=j; }
		    }
	    }
    }


    static final void find2LargestInRowStartingFrom ( IntegerMatrix matrix, int row, int startPoint, int[] pos, int[] val ) {
	pos[0]=-1;pos[1]=-1;val[0]=0;val[1]=0;

	matrix.checkRowIndex( row );
	matrix.checkColIndex( startPoint );

	int numOfRows = matrix.getNumRows();
	int numOfCols = matrix.getNumCols();

	for(int j=startPoint; j<numOfCols; j++)
	    {
		int k3=matrix.re[row][j];
		if( k3 != 0 )
		    {
			if( Math.abs(k3) > Math.abs(val[0]) )
			    {
				val[1]=val[0];pos[1]=pos[0];val[0]=k3;pos[0]=j;
			    } else if( Math.abs(k3) > Math.abs(val[1]) ) {  val[1]=k3; pos[1]=j; }
		    }
	    }
    }

    public static final IntegerMatrix invert( IntegerMatrix matrix1 ) {
	IntegerMatrix matrix2 = new IntegerMatrix ( matrix1.getNumRows() );

	invert( matrix1, matrix2 );

	return matrix2;
    }


    public static final void invert ( IntegerMatrix matrix1, IntegerMatrix matrix2 ) {
	if( !matrix1.isSquared() )
	    throw new IllegalArgumentException("can not invert non squared matrix");

	int dimension = matrix1.getNumRows();

	IntegerMatrix A = new IntegerMatrix ( matrix1 );

	matrix2.assignId( dimension );

        int[] pos = new int[2];
        int[] val = new int[2];

	for(int i=0;i<dimension;i++)
	    {
		while(true)
		    {

			find2LargestInColStartingFrom ( A, i, i, pos, val );
			if( pos[0] == -1 )
			    {
				throw new IllegalArgumentException("matrix is singular");
			    } else
				{
				    if ( pos[1] == -1 && pos[0] == i) break;
				    if ( pos[1] == -1 && pos[0] != i)
					{
					    swapRows( A, i, pos[0] );
					    swapRows( matrix2, i, pos[0] );
					}
				    if ( pos[1] != -1)
					{
					    int k1=-val[0]/val[1];
					    addRow( A, pos[0], pos[1], k1 );
					    addRow( matrix2, pos[0], pos[1], k1 );
					}
				}
		    }
	    }

	int det=1;

	for(int i=0;i<dimension;i++)
	    det*=A.re[i][i];

	if( det<0 )
	    det *= -1;

	for(int i=1;i<dimension;i++)
	    {
		int k1=A.re[i][i];
		for(int j=0;j<i;j++)
		    {
			int k2=-A.re[j][i];
			multRow ( A  , j , k1 );
			multRow ( matrix2  , j , k1 );
			addRow( A, j, i, k2 );
			addRow( matrix2, j, i, k2 );
		    }
	    }
	for(int i=0;i<dimension;i++)
	    {
		int k1=det/A.re[i][i];
		multRow ( A        , i , k1 );
		multRow ( matrix2  , i , k1 );
	    }
    }

    /////****


    //*****

     static final void moduloTwoRes(IntegerMatrix mat1, IntegerMatrix mat2) {
       int numOfRows = mat1.getNumRows();
       int numOfCols = mat1.getNumCols();

       mat2.newSize(numOfRows, numOfCols);

       for (int i = 0; i < numOfRows; i++)
         for (int j = 0; j < numOfCols; j++) {
           int tmp = mat1.re[i][j];
           tmp -= 2 * (tmp / 2);
           if (tmp < 0)
             tmp *= -1;
           mat2.re[i][j] = tmp;
         }
     }

    static public final int rankModuloTwo( IntegerMatrix tau ) {

	int zero_space = 0;
	int genus      = tau.getNumRows();

	int[] pos = new int[2];
	int[] el  = new int[2];

	IntegerMatrix A = new IntegerMatrix ( tau );

	for(int i=0 ; i+zero_space < genus ; i++) {
	    while( true ) {

		moduloTwoRes ( A , A );

		find2LargestInColStartingFrom ( A, i, i, pos, el );

		if( pos[0] == -1 ) {
		    zero_space++;

		    if(i + zero_space >=  genus) break;

		    moveColAndRowToEnd ( A, i );

		} else {

		    if( A.get( i, i ) == 1) {

			if(pos[1] == -1) break;

			for( int k=i+1; k< genus; k++) {

			    int l = A.get( k, i );

			    addColAndRow ( A, k, i , -l);
			    moduloTwoRes ( A , A  );
			}
			// A.print( " mod2 tau " );
			break;

		    } else {

			swapColsAndRows ( A, i+1, pos[0] );

			//  A.print( " mod2 tau ");

			for(int k=i+2; k< genus; k++) {

			    int l = A.get( k, i );

			    addColAndRow ( A, k, i+1 , -l );
			    moduloTwoRes ( A , A );
			}

			i++;

			for( int k=i+1; k< genus; k++) {

			    int l = A.get(k,i);
			    addColAndRow ( A, k, i-1 , -l);
			    moduloTwoRes ( A , A);
			}
			// A.print(" mod2 tau ");
			break;
		    }
		}
	    }
	}

	return  genus - zero_space;
    }

}









