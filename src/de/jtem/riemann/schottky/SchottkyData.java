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

package de.jtem.riemann.schottky;

import java.io.Serializable;

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;

/**
 * @version 0.5
 * @author schmies
 *
 */
public class SchottkyData
    implements Serializable {

  private static final long serialVersionUID = 1L;

  protected SchottkyGroupElement [] generator;

  protected SchottkyGroupElement [] generatorInv;

  /**
   ** Number of generators of schottky group.
   **/
  final int numGenerators;

  /**
   * Radius of ismetric circles.
   */
  double[] radius;

  /**
   * Centers of isometric cirlces.
   */
  Complex[][] center;

  /**
   * Fixpoints of loxodromic generators.
   */
  Complex[][] fixpoint;

  /**
   * Field.Complex factor of loxodromic generators.
   */
  Complex[] mu;


  double[] uniformizationData;

  public SchottkyData() {
    this(2);
  }

  public SchottkyData(int numOfGenerators) {
    this(getDefaultUniformizationData(numOfGenerators));
  }

  public SchottkyData(double[] uniformizationData) {

    if (uniformizationData.length % 6 != 0) {
      throw new IllegalArgumentException(
          " uniformizatinoData has wrong length ");
    }

    numGenerators = uniformizationData.length / 6;

    init();

    System.arraycopy(uniformizationData, 0,
                     this.uniformizationData, 0, uniformizationData.length);

    updateFixPointsAndMusFromDataArray();
    updateGeneratorsFromFixPointsAndMus();
    updateCircles();
  }

  private void init() {

    generator = new SchottkyGroupElement[numGenerators];
    generatorInv = new SchottkyGroupElement[numGenerators];

    for ( char i = 0; i < numGenerators; i++) {
      generator[i] = new SchottkyGroupElement();
      generatorInv[i] = new SchottkyGroupElement();

      generator[i].left = generator[i].right = i;
      generator[i].leftIsInvert = generator[i].rightIsInvert = 0;
      generator[i].wordLength = 1;

      generatorInv[i].left = generatorInv[i].right = i;
      generatorInv[i].leftIsInvert = generatorInv[i].rightIsInvert = 1;
      generatorInv[i].wordLength = 1;
    }

    radius = new double[numGenerators];
    uniformizationData = new double[numGenerators * 6];

    center = new Complex[numGenerators][2];
    fixpoint = new Complex[numGenerators][2];

    mu = new Complex[numGenerators];

    for (int i = 0; i < numGenerators; i++) {

      center[i][0] = new Complex();
      center[i][1] = new Complex();

      fixpoint[i][0] = new Complex();
      fixpoint[i][1] = new Complex();

      mu[i] = new Complex();
    }
  }

  public double[] getUniformizationData() {
    return (double[]) uniformizationData.clone();
  }

  public Moebius getGenerator(int anIndex) {
    return new Moebius(generator[anIndex]);
  }

  public double[] getRadius() {
    return radius;
  }

  public double getRadius(int anIndex) {
    return radius[anIndex];
  }

  public Complex[][] getCenters() {
    return center;
  }

  public Complex getCenterOfCircle(int anIndex, boolean isPrime) {
    return center[anIndex][isPrime ? 1 : 0];
  }

  public Complex getCenterOfCircle(int anIndex) {
    return center[anIndex][0];
  }

  public void getCenterOfCircle(int anIndex, boolean isPrime, Complex center ) {
  	center.assign( this.center[anIndex][isPrime ? 1 : 0] );
  }

  public void getCenterOfCircle(int anIndex, Complex center ) {
  	center.assign( this.center[anIndex][0] );
  }
  
  
  final IllegalArgumentException generatorsIntersect =
      new IllegalArgumentException("Generators intersect");

  static double[] getDefaultUniformizationData(int genus) {

    double[] data = new double[6 * genus];

    double muRe = 0.01;

    for (int i = 0, j = 0; i < genus; i++, muRe /= 5) {

      data[j++] = i + 1;
      data[j++] = 0;
      data[j++] = -i - 1;
      data[j++] = 0;
      data[j++] = muRe;
      data[j++] = 0;
    }

    return data;
  }

  public int getNumGenerators() {
    return numGenerators;
  }

  public Complex getA(int anIndex) {
    return new Complex(fixpoint[anIndex][0]);
  }

  public void setA(int anIndex, Complex aComplex) {
    fixpoint[anIndex][0].assign(aComplex);
    updateFromFixPointsAndMus();
  }

  public Complex getB(int anIndex) {
    return new Complex(fixpoint[anIndex][1]);
  }

  public void setB(int anIndex, Complex aComplex) {
    fixpoint[anIndex][1].assign(aComplex);
    updateFromFixPointsAndMus();
  }

  public void getMu(int anIndex, Complex mu) {
    mu.assign(this.mu[anIndex]);
  }

  public Complex getMu(int anIndex) {
    return new Complex(mu[anIndex]);
  }

  public void setMu(int anIndex, Complex aComplex) {
    mu[anIndex].assign(aComplex);
    updateFromFixPointsAndMus(); ;
  }

  public double[] getDoubleArrayValue() {
    double[] data = new double[6 * numGenerators];

    getUniformizationData(data, 0);

    return data;
  }

  public int getDoubleArrayValueLength() {
    return 6 * numGenerators;
  }

  public void getValue(double[] data) {
    getUniformizationData(data, 0);
  }

  public void getUniformizationData(double[] data, int offset) {
    System.arraycopy(uniformizationData, 0, data, offset,
                     uniformizationData.length);
  }

  public void setUniformizationData(double[] data) {
    setDoubleArrayParameter(data, 0);
  }

  final private boolean equalToUniformizationData( double[] data, int offset ) {
	for( int i=0, j=offset; i<uniformizationData.length; i++, j++ )
		  if( uniformizationData[i]!=data[j])
		  	return false;
  		return true;
  }
  
  public void setDoubleArrayParameter(double[] data, int offset) {
  	
  	if( equalToUniformizationData( data, offset ) )
  		return;
  		
    System.arraycopy(data, offset, uniformizationData, 0,
                     uniformizationData.length);

    updateFixPointsAndMusFromDataArray();
    updateGeneratorsFromFixPointsAndMus();
    update();
  }

  public void set(
      Complex[] A,
      Complex[] B,
      Complex[] Mu) {

    if (A.length != numGenerators
        || B.length != numGenerators
        || Mu.length != numGenerators) {
      throw new IllegalArgumentException("array(s) have wrong length");
    }

    for (int i = 0; i < numGenerators; i++) {
      fixpoint[i][0].assign(A[i]);
      fixpoint[i][1].assign(B[i]);
      mu[i].assign(Mu[i]);
    }

    updateFromFixPointsAndMus();
  }

  /**
   * @deprecated
   * @param data
   */
  public void set(double[] data) {
    setUniformizationData(data);
  }

  public void set(Moebius[] generator) {

    if (generator.length != numGenerators) {
      throw new IllegalArgumentException("array(s) have wrong length");
    }

    for (int i = 0; i < numGenerators; i++) {
      this.generator[i].assign(generator[i]);
      generatorInv[i].assignInvert(this.generator[i]);
    }

    updateFixPointsAndMusFromGenerators();
    updateDataArrayFromFixPointsAndMus();
    update();
  }

  void updateGeneratorsFromFixPointsAndMus() {

    for (int i = 0; i < numGenerators; i++) {
      generator[i].assign(fixpoint[i][0], fixpoint[i][1], mu[i]);
      generatorInv[i].assignInvert(generator[i]);
      generatorInv[i].parent = generator[i].parent;
    }
  }

  void updateFromFixPointsAndMus() {
    updateDataArrayFromFixPointsAndMus();
    updateGeneratorsFromFixPointsAndMus();
    update();
  }

  void updateDataArrayFromFixPointsAndMus() {

    for (int i = 0, j = 0; i < numGenerators; i++) {
      uniformizationData[j++] = fixpoint[i][0].re;
      uniformizationData[j++] = fixpoint[i][0].im;
      uniformizationData[j++] = fixpoint[i][1].re;
      uniformizationData[j++] = fixpoint[i][1].im;
      uniformizationData[j++] = mu[i].re;
      uniformizationData[j++] = mu[i].im;
    }
  }

  void updateFixPointsAndMusFromDataArray() {

    for (int i = 0, j = 0; i < numGenerators; i++) {

      fixpoint[i][0].assign(uniformizationData[j++], uniformizationData[j++]);
      fixpoint[i][1].assign(uniformizationData[j++], uniformizationData[j++]);

      mu[i].assign(uniformizationData[j++], uniformizationData[j++]);
    }
  }

  private Complex[] eigenvalues = new Complex[2];

  void updateFixPointsAndMusFromGenerators() {

    for (int i = 0; i < numGenerators; i++) {
      generator[i].getEigenValues(eigenvalues);
      mu[i].assignSqr(eigenvalues[0]);
      generator[i].getFixPoints(fixpoint[i]);
    }
  }

  void updateCircles() {
    for (int i = 0; i < numGenerators; i++) {
      radius[i] = generator[i].getCircles(center[i]);
    }
  }

  /**
    * Querries if this generates a classical Schottky group.
    * @return true iff this generates a classical Schottky group.
    */
   public boolean isClassical() {

     for( int i=0, N=0; i<2; i++ ) {
       for( int n=0; n<numGenerators; n++, N++ ) {

         final Complex c1 = center[n][i];
         final double r1 = radius[n];

         for( int j=0, M=0; j<2; j++ ) {
           for (int m = 0; m < numGenerators && M<N; m++, M++) {

             final Complex c2 = center[m][j];
             final double r2 = radius[m];

             if( c1.distSqr( c2 ) < (r1+r2)*(r1+r2) )
               return false;
           }
         }
       }
     }

     return true;
   }

  void update() {
    updateCircles();
  }

  /** 
   * Checks whether a point <code>p</code> is inside the
   * fundamental domain <code>F</code>.
   * The boundary of <code>F</code>, i.e. the isometric circles
   * are considered to be part of <code>F</code>,
   * @see #isInFundamentalDomain(Complex, double)
   */
  public boolean isInFundamentalDomain(Complex p ) {
	  return isInFundamentalDomain( p, -1e-12 );
  }
  /**
   * Checks whether a point <code>p</code> is inside the 
   * fundamental domain <code>F</code>.
   * With relTol you provide a signed relative tolerance, which means
   * that the creteria scale the radii by a factor of (1-relTol).
   * Thus a positive tolerance will enlarge the fundamental domain by the
   * given relative tolerance.
   */
  public boolean isInFundamentalDomain(Complex p, double relTol) {

    for (int i = 0; i < numGenerators; i++) {

      final double thresh = radius[i] * radius[i] * ( 1 - 2*relTol ) ;

      for (int j = 0; j < 2; j++) {

        if (center[i][j].distSqr(p) < thresh )
          return false;
      }
    }

    return true;
  }

  /** Returns distance of P to boundary of fundamental domain F.
   * Thus the result is only zero if P lays in one
   * of the 2N isomorphics circles which forms the boundary of F. */
  public double distToBoundaryOfFundamentalDomain(Complex P) {

    double minDist = Double.MAX_VALUE;

    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < numGenerators; i++) {

        double aDist = radius[i] - center[i][j].dist(P);

        if (aDist > 0) {
          return aDist;
        }

        if ( -aDist < minDist) {
          minDist = -aDist;
        }
      }
    }

    return minDist;
  }

}