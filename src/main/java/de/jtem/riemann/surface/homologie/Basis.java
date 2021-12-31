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
import de.jtem.blas.IntegerVector;

/** 
    All members of this class are defined during instantiation and stay unchanged afterwards:</br>

    <i>numOfCycles</i></br>
    <i>statusOfEdges</i></br>
    <i>dependenceMatrix</i></br>
    <i>firstCycleBasis</i></br>
    <i>firstCycleIntersections</i></br>

  */

public class Basis extends Transform implements Serializable, Cloneable {

    private static final long serialVersionUID = 1L;

    int numOfCycles;

    IntegerMatrix dependenceMatrix;
    IntegerMatrix firstCycleBasis;
    IntegerMatrix firstCycleIntersections;

    IntegerVector statusOfEdges;

 
    public void assign( Basis aTransform ) {
	super.assign( aTransform );

	numOfCycles = aTransform.numOfCycles;

	if( aTransform.dependenceMatrix == null ) dependenceMatrix = null;
	else if( dependenceMatrix == null ) 
	    dependenceMatrix = new IntegerMatrix ( aTransform.dependenceMatrix  );
	else dependenceMatrix.assign(  aTransform.dependenceMatrix  );

	if( aTransform.firstCycleBasis == null ) firstCycleBasis = null;
	else if( firstCycleBasis == null )  
	    firstCycleBasis = new IntegerMatrix ( aTransform.firstCycleBasis  );
	else firstCycleBasis.assign(  aTransform.firstCycleBasis  );

	if( aTransform.firstCycleIntersections == null ) firstCycleIntersections = null;
	else if( firstCycleIntersections == null ) 
	    firstCycleIntersections = new IntegerMatrix ( aTransform.firstCycleIntersections  );
	else firstCycleIntersections.assign(  aTransform.firstCycleIntersections  );

	if( aTransform.statusOfEdges == null ) statusOfEdges = null;
	else if( statusOfEdges == null ) 
	    statusOfEdges = new IntegerVector ( aTransform.statusOfEdges  );
	else statusOfEdges.assign(  aTransform.statusOfEdges  );
    }

    public Object clone() {
	return new Basis( this );
    }

    /**
       * Get the value of numOfCycles.
       * @return Value of numOfCycles.
       */
    public int getNumOfCycles() {return numOfCycles;}
    
    /**
       * Get the value of dependenceMatrix.
       * @return Value of dependenceMatrix.
       */
    public IntegerMatrix getDependenceMatrix() {return dependenceMatrix;}
    
    /**
       * Get the value of firstCycleBasis.
       * @return Value of firstCycleBasis.
       */
    public IntegerMatrix getFirstCycleBasis() {return firstCycleBasis;}

    /**
       * Get the value of firstCycleIntersections.
       * @return Value of firstCycleIntersections.
       */
    public IntegerMatrix getFirstCycleIntersections() {return firstCycleIntersections;}

    /**
       * Get the value of statusOfEdges.
       * @return Value of statusOfEdges.
       */
    public IntegerVector getStatusOfEdges() {return statusOfEdges;}
    
    protected Basis() {
    }

    public Basis( Basis aBasis ) {       
	assign( aBasis );
    }

    public Basis( IntegerMatrix monodrom, IntegerMatrix singular, IntegerVector disting ) {

	super( monodrom, singular, disting );

	int i, j, k, k1, k2=-1, k3=-1, k4=-1, k5, k6, k7, k8;
	int current_edge, current_g, current_cell, numOfMarkedPoints;
	int s_p1, s_p2;
	int numIntPoints, intPoint1, intPoint2;
	int in_v1, in_v2, out_v1, out_v2, inter_num;
	int thisLevelFound, level;

	IntegerVector gtmp                 = new IntegerVector ( numOfSheets );
	IntegerVector dependenceVector     = new IntegerVector ( numOfEdges );
	IntegerMatrix boundary2CellsMatrix = new IntegerMatrix ( numOfCells, numOfEdges );

	dependenceMatrix     = IntegerMatrix.id( numOfEdges );
 
	/* Now we calculate the relations generated by 2-cells */

	current_cell=-1;

	for( i=0; i<numOfBranchPoints; i++ ) 
	    {
		gtmp.assign( 0 );
		
		for( j=0; j<numOfSheets; j++ )
		    {
			current_g=g.re[i][j];
			if((current_g!=j)&&( gtmp.re[j]==0)) 
			    {	
				current_cell++;
				current_g=j;
				while(true)
				    {
					gtmp.re[current_g]=1;
					current_edge=gl.re[i][current_g];
					if(current_edge>-1) boundary2CellsMatrix.re[current_cell][current_edge]=1; 
					current_g=g.re[i][current_g];
					if(current_g==j) break;
				    }   
			    }
		    }
	    }

	gtmp.assign( 0 );


	for( j=0; j<numOfSheets; j++ )
	    {  
		if(gtmp.re[j]==0)
		    {
			current_cell++;
			current_g=j;	
			while(true)
			    {
				gtmp.re[current_g]=1;
				k=current_g;
				for( i=0; i<numOfBranchPoints; i++ )
				    {
					current_edge=gl.re[i][k];
					if(current_edge>-1) boundary2CellsMatrix.re[current_cell][current_edge]=1; 
					k=g.re[i][k];
				    }
				current_g=k;
				if(current_g==j) break;
			    }
		    }      
	    }
  

	for( i=0; i< numOfCells; i++) 
	    {
		current_edge=-1;
		for ( j=numOfEdges-1; j>=0; j--)
		    {
			if(boundary2CellsMatrix.re[i][j]!=0)
			    {
				dependenceVector.re[j]=1;
				for( k=0; k< numOfEdges; k++)
				    {
					k1=dependenceMatrix.re[k][j]/boundary2CellsMatrix.re[i][j];
					for( k2=0; k2< numOfEdges; k2++) 
					    dependenceMatrix.re[k][k2]-=k1*boundary2CellsMatrix.re[i][k2];
				    }  
				for( k=i+1; k< numOfCells; k++)
				    {
					k1=boundary2CellsMatrix.re[k][j]/boundary2CellsMatrix.re[i][j];
					for( k2=0; k2< numOfEdges; k2++) 
					    boundary2CellsMatrix .re[k][k2]-=k1*boundary2CellsMatrix.re[i][k2];
				    }
				break;
			    }
		    }
	    }

	if( debug )
	    dependenceMatrix.print( "Dependence  Matrix" );

	/* Now we calculte the spanning tree inside the graph of edges */

	IntegerVector levelVector       = new IntegerVector ( numOfSheets );
	IntegerVector prevLevelVector   = new IntegerVector ( numOfSheets );
	IntegerVector edgeDownNumber    = new IntegerVector ( numOfSheets );
	IntegerVector edgeOrientation   = new IntegerVector ( numOfSheets );
	
	statusOfEdges     = new IntegerVector ( numOfPaths );

	for ( j=0; j <numOfSheets; j++)
	    { 
		levelVector.re[j]=-1;
		prevLevelVector.re[j]=-2;
		edgeOrientation.re[j]=0;
	    }
	for ( j=0; j <numOf_g_Edges; j++)
	    { 
		k=dependenceVector.re[j];
		statusOfEdges.re[j]=-k;
	    }
	
	level=0;
	levelVector.re[0]=level;
	prevLevelVector.re[0]=-1;  
	edgeDownNumber.re[0]=-1;
	numOfMarkedPoints=1;

	do
	    {
		thisLevelFound=0;
		for( k1=0; k1<numOf_g_Edges; k1++ ) 
		    {
			k2= statusOfEdges.re[k1];
			if(k2==0)
			    {
				i=edgeBranchPoint.re[k1];
				j=edgeStartPoint.re[k1];	      
				k=edgeEndPoint.re[k1];	      
				k3=levelVector.re[j];
				k4=levelVector.re[k];
				if(k3==level&&k4==-1)
				    {
					thisLevelFound++;
					levelVector.re[k]=level+1;
					prevLevelVector.re[k]=j;
					edgeDownNumber.re[k]=k1;
					statusOfEdges.re[k1]=1;
					edgeOrientation.re[k]=1;		  
					numOfMarkedPoints++;
				    }
				if(k4==level&&k3==-1)
				    {
					thisLevelFound++;
					levelVector.re[j]=level+1;
					prevLevelVector.re[j]=k; 
					edgeDownNumber.re[j]=k1;
					statusOfEdges.re[k1]=1;
					edgeOrientation.re[j]=-1;  
					numOfMarkedPoints++;
				    }
			    }
		    }
		level++;
	    } while( thisLevelFound!=0 );

	for( i=0; i<numOfSingularPoints; i++ ) 
	    {
		for ( j=0; j <numOfSheets; j++)
		    { 
			k=h.re[i][j];
			if(k>-1) statusOfEdges.re[k]=2;
		    }
	    }
	
	for( i=0; i<numOfDistinguishedPoints; i++ ) 
	    {
		k=i+numOf_g_Edges+numOf_h_Edges;
		statusOfEdges.re[k]=4;
	    } 

	/* The spanning tree is constructed */ 
	
	if( debug )
	    System.out.println( " numOfMarkedPoints = "+ numOfMarkedPoints + "    level= " + level );
	
	if( numOfMarkedPoints != numOfSheets )
	    throw new IllegalArgumentException( "Surface is not connected" );

	if( debug ){
	    levelVector.print( "Level vector:" );
	    prevLevelVector.print( "Previous level vector:" );
	    edgeDownNumber.print( "Edge down number:" );
	    edgeOrientation.print( "Edge orientation:" );
      

	//PrintVectorEnumeratedBy(edgeBranchPoint, edgeDownNumber, numOfSheets, 3 ,   " Edge branch point     =   " ) ;
	//PrintVectorEnumeratedBy(edgeStartPoint,  edgeDownNumber, numOfSheets, 3 ,   " Edge start  point     =   " ) ;
	//PrintVectorEnumeratedBy(edgeEndPoint,    edgeDownNumber, numOfSheets, 3 ,   " Edge end  point       =   " ) ;

	    statusOfEdges.print( "Status of Edges:" );
	}
	/* Now we calculate the  cycles associated with independent edges outside the spanning tree  */
	
	numOfCycles=0;
	for ( j=0; j <numOf_g_Edges; j++)
	    {
		k=statusOfEdges.re[j];
		if( k==0) numOfCycles++;
	    }

	if( debug )
	    System.out.println("number of cycles = " + numOfCycles );

	IntegerVector basicEdgeNumbers = new IntegerVector ( numOfCycles );

	current_g=-1;
	for ( j=0; j <numOfEdges; j++)
	    {
		k=statusOfEdges.re[j];
		if( k==0) 
		    {
			current_g++;
			basicEdgeNumbers.re[current_g]=j;
		    }
	    }

	if( debug )
	    basicEdgeNumbers.print("Basic Edges:" );

	IntegerMatrix basisPathVertices      = new IntegerMatrix ( numOfCycles, numOfSheets );
	IntegerMatrix basisPathEdges         = new IntegerMatrix ( numOfCycles, numOfSheets );
	IntegerMatrix basisPathOrientations  = new IntegerMatrix ( numOfCycles, numOfSheets );
	
	IntegerVector basisPathLength        = new IntegerVector ( numOfCycles );

	firstCycleBasis        = new IntegerMatrix ( numOfCycles, numOfEdges );

	for ( i=0; i <numOfCycles; i++)
	    {
		for ( j=0; j <numOfSheets; j++)
		    {
			basisPathVertices.re[i][j]=-1;
			basisPathEdges.re[i][j]=-1; 
			basisPathOrientations.re[i][j]=0; 
		    }
		for ( j=0; j <numOfEdges; j++) firstCycleBasis.re[i][j]=0; 
		k=basicEdgeNumbers.re[i];
		s_p1=0;s_p2=numOfSheets-1;

		k2=edgeStartPoint.re[k];
		k1=edgeEndPoint.re[k];
		k3=levelVector.re[k1];
		k4=levelVector.re[k2];
    
		basisPathVertices.re[i][s_p2]=k2;
		basisPathEdges.re[i][s_p2]=k; 
		basisPathOrientations.re[i][s_p2]=1;   
		firstCycleBasis.re[i][k]=1;
		current_g=1;
		s_p2--;
		do
		    {
			k3=levelVector.re[k1];
			k4=levelVector.re[k2];  
			if(k3>=k4)
			    {
				k=edgeDownNumber.re[k1];
				basisPathEdges.re[i][s_p1]=k;
				basisPathOrientations.re[i][s_p1]=-edgeOrientation.re[k1];
				firstCycleBasis.re[i][k]=basisPathOrientations.re[i][s_p1];
				basisPathVertices.re[i][s_p1]=k1;
				k1=prevLevelVector.re[k1];
				s_p1++;
				current_g++;
	
			    }
			if(k3<=k4)
			    {	
				k=edgeDownNumber.re[k2];
				basisPathEdges.re[i][s_p2]=k;
				basisPathOrientations.re[i][s_p2]=edgeOrientation.re[k2];
				firstCycleBasis.re[i][k]=basisPathOrientations.re[i][s_p2];   
				k2=prevLevelVector.re[k2];
				basisPathVertices.re[i][s_p2]=k2;
				s_p2--;
				current_g++;
			    }
		    } while (k1!=k2);
		basisPathLength.re[i]=current_g;
		k1=numOfSheets-current_g;
		for(; s_p1+k1<numOfSheets; s_p1++)
		    {
			s_p2=s_p1+k1;
			basisPathVertices.re[i][s_p1]=basisPathVertices.re[i][s_p2];
			basisPathEdges.re[i][s_p1]=basisPathEdges.re[i][s_p2];
			basisPathOrientations.re[i][s_p1]= basisPathOrientations.re[i][s_p2];

		    }

	    }

	/* The first basis of cycles was calculated */
	
	if( debug ) {
	    basisPathLength.print( "Basis path lengths:" );

	    System.out.println("Path Vertices:");
	
	    for( i=0; i< numOfCycles; i++) 
		{
		    StringBuffer aLine = new StringBuffer();

		    k1=basisPathLength.re[i];
		    for ( j=0; j <k1; j++)
			{ 
			    k=basisPathVertices.re[i][j];
			    aLine.append("  "+k );
			}
		    System.out.println( aLine );
		}
	    System.out.println(" Path Edges \n");
	    
	    for( i=0; i<numOfCycles; i++) 
		{      
		    StringBuffer aLine = new StringBuffer();
		    
		    k1=basisPathLength.re[i];
		    for ( j=0; j <k1; j++)
			{ 
			    k=basisPathEdges.re[i][j];
			    k=k*basisPathOrientations.re[i][j];
			    aLine.append("  "+k );
			}
		    System.out.println( aLine );
		}
	
	    firstCycleBasis.print( "First Cycle basis" );
	}

	/* The basis of independent edges is calculated !!!!!!!!!!!!!!!!!!!!!!! */
	
	
	/* Now we calculate the intersection matrix !!!!!!!!!!!!!!!!!!!!!!!!!!! */

	firstCycleIntersections = new IntegerMatrix ( numOfCycles,  numOfCycles );
	
	for(i=0; i< numOfCycles; i++)
	    {      
		k5=basisPathLength.re[i];
		for(j=0; j< numOfCycles; j++)
		    {

			firstCycleIntersections.re[i][j]=-10;   

			if(i==j)  firstCycleIntersections.re[i][j]=0; 
			else
			    {
				k6=basisPathLength.re[j];  
				numIntPoints=0;
				intPoint1=-1;
				intPoint2=-1;
				for(k1=0; k1<k5; k1++)
				    {
					k7=basisPathVertices.re[i][k1];
					for(k2=0; k2<k6; k2++)
					    {
						k8=basisPathVertices.re[j][k2];	  
						if(k8==k7) 
						    {
							numIntPoints++;
							intPoint1=k7;
							break;
						    }      
					    }
					if(intPoint1>-1) break;
				    }
				if(numIntPoints>0) 
				    {	    
					for(k3=k5-1; k3>=0; k3--)
					    {
						k7=basisPathVertices.re[i][k3];
						for(k4=k6-1; k4>=0; k4--)
						    {
							k8=basisPathVertices.re[j][k4];	  
							if(k8==k7) 
							    {	
								intPoint2=k7;
								break;
							    }  
						    }
						if(intPoint2>-1) break;
					    }
					
					if(intPoint2!=intPoint1) numIntPoints++;
				    }
				if(numIntPoints==0) firstCycleIntersections.re[i][j]=0;
				else
				    {
					inter_num=0;
					k=k1;
					k7=basisPathEdges.re[i][k];
					k8=edgeBranchPoint.re[k7];
					out_v1=2*k8;
					k7=basisPathOrientations.re[i][k];
					if(k7<0) out_v1++;
					if(k==0) k=basisPathLength.re[i]-1;
					else k--;
					k7=basisPathEdges.re[i][k];
					k8=edgeBranchPoint.re[k7];
					in_v1=2*k8;
					k7=basisPathOrientations.re[i][k];
					if(k7>0) in_v1++;
					k=k2;
					k7=basisPathEdges.re[j][k];
					k8=edgeBranchPoint.re[k7];
					out_v2=2*k8;
					k7=basisPathOrientations.re[j][k];
					if(k7<0) out_v2++;
					if(k==0) k=basisPathLength.re[j]-1;
					else k--;
					k7=basisPathEdges.re[j][k];
					k8=edgeBranchPoint.re[k7];
					in_v2=2*k8;
					k7=basisPathOrientations.re[j][k];
					if(k7>0) in_v2++;
	  
					if(out_v1>in_v1)
					    {
						if((in_v2>in_v1)&&(in_v2<out_v1)) inter_num++;
						if((in_v2<in_v1)||(in_v2>out_v1)) inter_num--;
						if((out_v2>in_v1)&&(out_v2<out_v1)) inter_num--;
						if((out_v2<in_v1)||(out_v2>out_v1)) inter_num++;	      
					    }
					
					if(out_v1<in_v1)
					    {
						if((in_v2>out_v1)&&(in_v2<in_v1)) inter_num--;
						if((in_v2<out_v1)||(in_v2>in_v1)) inter_num++;
						if((out_v2>out_v1)&&(out_v2<in_v1)) inter_num++;
						if((out_v2<out_v1)||(out_v2>in_v1)) inter_num--;	      
					    }


					//printf(" i = %3i j = %3i number of int %3i   ",i,j,numIntPoints);
					//printf(" in1 = %3i out1 = %3i in2 = %3i out2 = %3i   ",in_v1,out_v1,in_v2,out_v2);

					if(numIntPoints>1)
					    {
						
						k=k3;
						k7=basisPathEdges.re[i][k];
						k8=edgeBranchPoint.re[k7];
						out_v1=2*k8;
						k7=basisPathOrientations.re[i][k];
						if(k7<0) out_v1++;
						if(k==0) k=basisPathLength.re[i]-1;
						else k--;
						k7=basisPathEdges.re[i][k];
						k8=edgeBranchPoint.re[k7];
						in_v1=2*k8;
						k7=basisPathOrientations.re[i][k];
						if(k7>0) in_v1++;
						k=k4;
						k7=basisPathEdges.re[j][k];
						k8=edgeBranchPoint.re[k7];
						out_v2=2*k8;
						k7=basisPathOrientations.re[j][k];
						if(k7<0) out_v2++;
						if(k==0) k=basisPathLength.re[j]-1;
						else k--;
						k7=basisPathEdges.re[j][k];
						k8=edgeBranchPoint.re[k7];
						in_v2=2*k8;
						k7=basisPathOrientations.re[j][k];
						if(k7>0) in_v2++;
						if(out_v1>in_v1)
						    {
							if((in_v2>in_v1)&&(in_v2<out_v1)) inter_num++;
							if((in_v2<in_v1)||(in_v2>out_v1)) inter_num--;
							if((out_v2>in_v1)&&(out_v2<out_v1)) inter_num--;
							if((out_v2<in_v1)||(out_v2>out_v1)) inter_num++;	      
						    }
						
						if(out_v1<in_v1)
						    {
							if((in_v2>out_v1)&&(in_v2<in_v1)) inter_num--;
							if((in_v2<out_v1)||(in_v2>in_v1)) inter_num++;
							if((out_v2>out_v1)&&(out_v2<in_v1)) inter_num++;
							if((out_v2<out_v1)||(out_v2>in_v1)) inter_num--;	      
						    }
						
						//printf(" in1 = %3i out1 = %3i in2 = %3i out2 = %3i ",in_v1,out_v1,in_v2,out_v2);
					    }
					firstCycleIntersections.re[i][j]=inter_num/2;
				    }
				//printf("i = %3i j = %3i intersection number = %3i \n",i,j,inter_num);    */
			    }
		    }
	    }
	/*  The intersection matrix is found !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  */
    }

}










