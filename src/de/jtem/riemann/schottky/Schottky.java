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

import de.jtem.blas.ComplexMatrix;
import de.jtem.blas.ComplexVector;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import de.jtem.mfc.matrix.Complex2By2;

public class Schottky
    extends SchottkyData
    implements Serializable {

  private static final long serialVersionUID = 1L;

  private static double DEFAULT_ACCURACY;
  
  private static int initCapacity = 1000;

  int updateID = 0;

  double acc;

  double theta1, q1;

  final Complex zOfK2
      = new Complex(Double.NaN);

  double[][] k2;

  int numOfElements;

  /* value for series evaluablity test */
  int l = 2;
  
  /* value for series evaluablity test */
  double C = 0.75;
  
  /**
   * Throwing an RuntimeException when computation needs more
   * than this number of elements.
   */
  long maxNumOfElements = 200000;
  
  /**
   * Longest element of global used subset S.
   */
  int maxWordLength = 1;

  /**
   * innerCenter[n][i] and give the centers of the images of the 2N-1 circles
   * which are mapped by the generator[n][i] into the interior of the
   * isometric circle with center[n][(i+1)%2] and radius[n][i], which is the
   * the image of its counterpart with center[n][i].*/
  Complex[][][] innerCenter;

  /**
       * innerRadius[n] give the 2N-1 radii of the circles with innerCenter[n][0,1].
   */
  double[][][] innerRadius;

  double[][] distToCenter;

  final Complex A = new Complex();
  final Complex B = new Complex();
  final Complex H = new Complex();
  final Complex h = new Complex();
  final Complex tmp = new Complex();
  final Complex tmp2 = new Complex();
  final Complex sum = new Complex();
  final Complex diff = new Complex();

  boolean useFancyError = true;

  AbelianDifferential abelianDifferential;

  AbelianIntegral abelianIntegral;

  PeriodMatrix periodMatrix;

  Sigma sigma;

  public Schottky() {
    this(2);
  }

  public Schottky(int genus) {
    this(getDefaultUniformizationData(genus), DEFAULT_ACCURACY);
  }

  public Schottky( SchottkyData data ) {
	  this( data, DEFAULT_ACCURACY );
  }
  
  public Schottky( SchottkyData data, double accuracy ) {
	  this( data.getUniformizationData(), accuracy );
  }
  
  public Schottky(double[] uniformizationData) {
    this(uniformizationData, 1e-7);
  }

  public Schottky(double[] uniformizationData, double accuracy) {

    super(uniformizationData);

    init();

    this.acc = accuracy;

    update();
  }

  private void init() {

    innerCenter = new Complex[numGenerators][2][2 * numGenerators - 1];
    innerRadius = new double[numGenerators][2][2 * numGenerators - 1];

    for (int i = 0; i < numGenerators; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2 * numGenerators - 1; k++) {
          innerCenter[i][j][k] = new Complex();
        }
      }
    }
    
    distToCenter = new double[numGenerators][2];

    k2 = new double[2][numGenerators];

    id = new SchottkyGroupElement();

    id.child = new SchottkyGroupElement[2 * numGenerators];

    id.norm = Double.POSITIVE_INFINITY;

    for (int i = 0, j = 0; i < numGenerators; i++) {
      generator[i].parent = generatorInv[i].parent = id;

      id.child[j++] = generator[i];
      id.child[j++] = generatorInv[i];

      generator[i].norm = Double.POSITIVE_INFINITY;
      generatorInv[i].norm = Double.POSITIVE_INFINITY;
    }

    computeNumOfElementsWithWordLength();
    computeNumOfElementsOfCosetWithWordLength();

    abelianDifferential = new AbelianDifferential(this);
    abelianIntegral = new AbelianIntegral(this);
    periodMatrix = new PeriodMatrix(this);
    sigma = new Sigma(this);

  }

  void update() {

    super.update();

    if (numGenerators > 1 && !isClassical()) {
      throw generatorsIntersect;
    }

    updateID++;

    updateInnerCircles();

    zOfK2.assign(Double.NaN);

    theta1 = theta1();

    q1 = theta1 * theta1 * (2 * numGenerators - 1);

    if (q1 >= 1) {

      theta1 = theta2();

      q1 = theta1 * theta1 * (2 * numGenerators - 1);

      if (q1 >= 1) {
        throw new RuntimeException("cannot guarantee conversion");
      }
    }

    maxWordLength = 0;

    updateElementTree();

    abelianIntegral.update();
    periodMatrix.update();
    sigma.update();
  }

  public int getNumElements() {
    return numOfElements;
  }

  void updateInnerCircles() {

    for (int n = 0; n < numGenerators; n++) {
      for (int i = 0; i < 2; i++) {

        // sigma maps F into Circle with center[n][i]

        Moebius sigma = i == 0 ? generatorInv[n] : generator[n];

        for (int m = 0, k = 0; m < numGenerators; m++) {
          for (int j = 0; j < 2; j++) {

            if (m != n || i == j) {

              innerRadius[n][i][k] =
                  sigma.getRadiusOfMappedCircle(center[m][j], radius[m],
                                                innerCenter[n][i][k]);

              // check that it is in the target circle
              if (center[n][i].dist(innerCenter[n][i][k]) + innerRadius[n][i][k]
                  > radius[n]) {
                System.out.println("computeInnerCircles: fail target");
              }
              k++;
            }
          }
        }
      }
    }
  }

  public int getL() {
	  return l;
  }
  
  public void setL( int l ) {
	  if( l<1 ) 
		  throw new IllegalArgumentException( "value has to be positive!");
	  this.l = l;
  }
  
  public double getC() {
	  return C;
  }
  
  public void setC( double C ) {
	  if( C<=0 || C>=1 )
		  throw new IllegalArgumentException( "value has to be in ]0,1[!");
	  this.C = C;
  }
  
  public boolean isIntegralSeriesEvaluable() {
	  double theta = theta(l);
	  return theta * theta * (2 * numGenerators - 1) < C;
  }
  
  public boolean isDifferentialSeriesEvaluable() {
	  double kappa = kappa(l);
	  return (2 * numGenerators - 1) / kappa / kappa < C;
  }
  
  public boolean isSeriesEvaluable() {
	  return isIntegralSeriesEvaluable() || isDifferentialSeriesEvaluable();
  }
  
  public double getAccuracy() {
    return acc;
  }

  public void setAccuracy(double accuracy) {
  	if(acc == accuracy )
  		return;
    acc = accuracy;
    update();
  }

  void setAccuracy(int accuracy) {
    setAccuracy(Math.pow(10.0, - (double) accuracy));
  }

  public int getMaxWordLength() {
    return maxWordLength;
  }

  long[] numOfElementsOfCosetWithWordLength;

  final void computeNumOfElementsOfCosetWithWordLength() {
    final int length;

    if (numGenerators == 1) {
      length = 101; // hackish: just a very long word;
    }
    else {
      // determines maximal array length such that long integer
      // does not overflow
      length = 2 + (int) (
          Math.log(Long.MAX_VALUE / (2. * numGenerators - 2)) /
          Math.log(2. * numGenerators - 1));
    }

    numOfElementsOfCosetWithWordLength = new long[length];

    numOfElementsOfCosetWithWordLength[0] = 1;
    numOfElementsOfCosetWithWordLength[1] = 2 * numGenerators - 2;

    for (int i = 2; i < length; i++) {
      numOfElementsOfCosetWithWordLength[i] =
          numOfElementsOfCosetWithWordLength[i - 1] * (2 * numGenerators - 1);
    }
  }

  /**
   * Returns number of element of cosets G/G_n (G_n\G) with
   * given word length.
   */
  final public long getNumOfElementsOfCosetWithWordLength(final int length) {
    return numOfElementsOfCosetWithWordLength[length];
  }

  long[] numOfElementsWithWordLength;

  final void computeNumOfElementsWithWordLength() {
    final int length;

    if (numGenerators == 1) {
      length = 101; // hackish: just a very long word;
    }
    else {
      // determines maximal array length such that long integer
      // does not overflow
      length = 2 + (int) (
          Math.log(Long.MAX_VALUE / (2. * numGenerators)) /
          Math.log(2. * numGenerators - 1));
    }

    numOfElementsWithWordLength = new long[length];

    numOfElementsWithWordLength[0] = 1;
    numOfElementsWithWordLength[1] = 2 * numGenerators;

    for (int i = 2; i < length; i++) {
      numOfElementsWithWordLength[i] =
          numOfElementsWithWordLength[i - 1] * (2 * numGenerators - 1);
    }
  }

  /**
   * Returns number of element of schottky group G with
   * given word length.
   */
  final public long getNumOfElementsWithWordLength(final int length) {
    return numOfElementsWithWordLength[length];
  }

  SchottkyGroupElement id;

  final void updateElementTree() {

    id.updateID = updateID;
    computeImageOfFixpoints(id);

    numOfElements = 1;

    for (int i = 0; i < numGenerators; i++) {
      computeConstants(generator[i]);
      computeConstants(generatorInv[i]);
    }

  }

  final void updateElement(final SchottkyGroupElement element) {

    if (element.updateID == updateID) {
      return;
    }

    if (element.leftIsInvert == 1) {
      element.assignTimes(generatorInv[element.left], element.parent);
    }
    else {
      element.assignTimes(generator[element.left], element.parent);
    }

    computeConstants(element);
  }

  final void computeConstants(final SchottkyGroupElement element) {

    numOfElements++;

    if( numOfElements > maxNumOfElements ) {
    	System.out.println( "stoped computations with more than " + maxNumOfElements +
    			" number of group elements.");
    	throw new RuntimeException("too many elements");
    }
    
    element.updateID = updateID;
    element.getC(tmp);

    final double cAbsSqr = tmp.absSqr();

    element.norm = 1 / cAbsSqr;

    computeLOfInverse(element);
    computeImageOfFixpoints(element);

    if (maxWordLength < element.wordLength) {
      maxWordLength = element.wordLength;
    }
  }

  final Complex targetCenter = new Complex();

  final void computeLOfInverse( final SchottkyGroupElement element ) {

    if( element.LOfInverse == null ||
        element.LOfInverse.length != numGenerators ) {
      element.LOfInverse = new double[numGenerators];
    }

    final double radius = this.targetRadiusForInverseTransformation( element, targetCenter );

    for( int n=0; n<numGenerators; n++ ) {
      final Complex A = fixpoint[n][0];
      final Complex B = fixpoint[n][1];
      element.LOfInverse[n] = A.dist(B)
          / ( A.dist( targetCenter ) - radius )
          / ( B.dist( targetCenter ) - radius );
    }
  }

  final void computeImageOfFixpoints(final SchottkyGroupElement element) {

    if (element.imageOfA == null ||
        element.imageOfA.length != numGenerators) {
      element.imageOfA = new Complex[numGenerators];
      element.imageOfB = new Complex[numGenerators];

      for (int i = 0; i < numGenerators; i++) {
        element.imageOfA[i] = new Complex();
        element.imageOfB[i] = new Complex();
      }
    }

    for (int i = 0; i < numGenerators; i++) {
      element.applyTo(fixpoint[i][0], element.imageOfA[i]);
      element.applyTo(fixpoint[i][1], element.imageOfB[i]);
    }
  }

  final void createLeftChilds(SchottkyGroupElement element) {
    element.createLeftChilds(numGenerators);
  }

  /**
   * Prepares k2 for value z.
   * @param z
   */
  final double[][] k2(Complex z) {

    if (z.equals(zOfK2)) {
      return k2;
    }

    zOfK2.assign(z);

    for (int i = 0; i < 2; i++) {
      for (int n = 0; n < numGenerators; n++) {
        k2[i][n] = k2(i, n, z);
      }
    }
    return k2;
  }


  /**
   * Prepares k2 for value z.
   * @param z
   */
  final double[][] kIndexed(Complex z, int wordLength ) {

    if (z.equals(zOfK2)) {
      return k2;
    }

    zOfK2.assign(z);

    for (int i = 0; i < 2; i++) {
      for (int n = 0; n < numGenerators; n++) {
        k2[i][n] = Double.MAX_VALUE;
      }
    }

    if (updateID != updateID) {
      update();
    }

    eval_k( id, z, k2, wordLength );

    return k2;
  }


  final private void eval_k(SchottkyGroupElement element,
                              Complex z,
                              double [][] k, int wordLength) {

    if (element.wordLength > wordLength) {
      return;
    }

    if (element.updateID != updateID) {
      updateElement(element);
    }

    if (element.wordLength == wordLength) {
      k[element.leftIsInvert][element.left]
          = Math.min(k(element, z), k[element.leftIsInvert][element.left] );
    }

    if (element.child == null) {
      createLeftChilds(element);
    }

    final SchottkyGroupElement[] child = element.child;

    final int numOfChilds = child.length;

    for (int i = 0; i < numOfChilds; i++) {
      eval_k(child[i], z, k, wordLength);
    }
  }

  final private double eval_k(SchottkyGroupElement element,
                              Complex z,
                              double k, int wordLength) {

    if (element.wordLength > wordLength) {
      return k;
    }

    if (element.updateID != updateID) {
      updateElement(element);
    }

    if (element.wordLength == wordLength) {
      k = Math.min(k(element, z), k);
    }

    if (element.child == null) {
      createLeftChilds(element);
    }

    final SchottkyGroupElement[] child = element.child;

    final int numOfChilds = child.length;

    for (int i = 0; i < numOfChilds; i++) {
      k = eval_k(child[i], z, k, wordLength);
    }
    return k;
  }

  /**
   * Computes k for given z and word lengths.
   */
  public final double k(Complex z, int wordLength) {

    if (updateID != updateID) {
      update();
    }

    return eval_k(id, z, Double.MAX_VALUE, wordLength);
  }


  public final void abelianDifferentialOf1stKind(Complex r,
                                                 Complex z, int n) {
    abelianDifferential.of1stKind(r, z, n, acc);
  }

  public final void abelianDifferentialOf1stKind(Complex r,
                                                 Complex z, int n,
                                                 double accuracy) {
    abelianDifferential.of1stKind(r, z, n, accuracy);
  }

  public final Complex abelianDifferentialOf1stKind(Complex z, int n) {
    Complex r = new Complex();
    abelianDifferential.of1stKind(r, z, n, acc);
    return r;
  }

  public final Complex abelianDifferentialOf1stKind(Complex z, int n,
      double accuracy) {
    Complex r = new Complex();
    abelianDifferential.of1stKind(r, z, n, accuracy);
    return r;
  }

  public final Complex abelianIntegralOf1stKind(Complex z, int n) {
    final Complex r = new Complex();
    abelianIntegralOf1stKind(r, z, n);
    return r;
  }

  public final Complex abelianIntegralOf1stKind(Complex z, int n, double acc) {
    final Complex r = new Complex();
    abelianIntegralOf1stKind(r, z, n, acc);
    return r;
  }

  public final void abelianIntegralOf1stKind(Complex r,
                                             Complex z, int n) {
    abelianIntegral.eval(r, z, n, acc);
  }

  // methods depending on m,n use tmp,tmp2,sum
  public final void abelianIntegralOf1stKind
      (final Complex r,
       final Complex z,
       final int n, final double accuracy) {
    abelianIntegral.eval(r, z, n, accuracy);
  }

  /** gamma is differential of tird kind having simple poles
      at A and B of residue -1 and  1 respectively
   **/
  public final void abelianDifferentialOf3rdKind(Complex r, Complex z,
                                                 Complex A, Complex B) {
    abelianDifferentialOf3rdKind(r, z, A, B, acc);
  }

  // Has simple poles at A and B with residues -1 and 1
  public final void abelianDifferentialOf3rdKind(Complex r, Complex z,
                                                 Complex A, Complex B,
                                                 double accuracy) {
    abelianDifferential.of3rdKind(r, z, A, B, accuracy);
  }

  /** gamma is differential of tird kind having simple poles
      at A and B of residue -1 and  1 respectively
   **/
  public final Complex abelianDifferentialOf3rdKind(Complex z,
      Complex A, Complex B) {
    Complex r = new Complex();
    abelianDifferentialOf3rdKind(r, z, A, B, acc);
    return r;
  }

  // Has simple poles at A and B with residues -1 and 1
  public final Complex abelianDifferentialOf3rdKind(Complex z,
      Complex A, Complex B,
      double accuracy) {
    Complex r = new Complex();
    abelianDifferential.of3rdKind(r, z, A, B, accuracy);
    return r;
  }

  public final Complex abelianIntegralOf3rdKind(Complex z,
                                                Complex A, Complex B) {
    Complex r = new Complex();
    abelianIntegral.of3rdKind(r, z, A, B, acc);
    return r;
  }

  public final Complex abelianIntegralOf3rdKind(Complex z,
                                                Complex A, Complex B,
                                                double acc) {
    Complex r = new Complex();
    abelianIntegral.of3rdKind(r, z, A, B, acc);
    return r;
  }

  public final void abelianIntegralOf3rdKind(Complex r,
                                             Complex z,
                                             Complex A,
                                             Complex B) {
    abelianIntegral.of3rdKind(r, z, A, B, acc);
  }

  public final void abelianIntegralOf3rdKind(Complex r,
                                             Complex z,
                                             Complex A,
                                             Complex B, double acc) {
    abelianIntegral.of3rdKind(r, z, A, B, acc);
  }

  
  final public void getV(final ComplexVector V, final double accuracy) {
    V.newSize(numGenerators);
    for (int n = 0; n < numGenerators; n++) {
      sigma.eval(tmp, n, accuracy);
      V.set(n, tmp);
    }
  }

  final public void getV(final ComplexVector V) {
    getV(V, acc);
  }

  final public ComplexVector getV(final double accuracy) {
    ComplexVector V = new ComplexVector(numGenerators);
    getV(V, accuracy);
    return V;
  }

  final public ComplexVector getV() {
    return getV(acc);
  }

  final void V(Complex r, int n) {
    sigma.eval(r, n, acc);
  }

  public final Complex V(int n, double accuracy) {
    Complex r = new Complex();
    sigma.eval(r, n, accuracy);
    return r;
  }


  final public void getV( int k, final ComplexVector V, final double accuracy) {
  	if( k<1) throw new IllegalArgumentException( "only positive powers are permitted"); 
  	V.newSize(numGenerators);
  	for (int n = 0; n < numGenerators; n++) {
  		sigma.evalPow(tmp, n, k, accuracy);
  		V.set(n, tmp);
  	}
  }
  

  final public ComplexVector getV( final int k, final double accuracy) {
  	ComplexVector V = new ComplexVector(numGenerators);
  	getV(k, V, accuracy);
  	return V;
  }

  final public ComplexVector getV( int k ) {
  	return getV( k, acc);
  }

  final void V(Complex r, int k, int n) {
  	sigma.evalPow(r, k, n, acc);
  }

  public final Complex V(int k, int n, double accuracy) {
  	Complex r = new Complex();
  	sigma.evalPow(r, n, k, accuracy);
  	return r;
  }
  
  
  public void sigma(Complex r, Complex z, Complex w, double acc) {
    sigma.eval(r, z, w, acc);
  }

  /**
   * sigma.
   * @param r result with default accuracy on output
   * @param z postion of evaluation
   * @param w postion of evaluation
   */
  public void sigma(Complex r, Complex z, Complex w) {
    sigma.eval(r, z, w, acc);
  }

  /**
   * sigma with power of k.
   * @param r result on output
   * @param z postion of evaluation
   * @param k power
   * @param acc accuracy of result
   */
  public void sigma(Complex r, Complex z, Complex w, int k, double acc) {
    sigma.eval(r, z, w, k, acc);
  }

  /**
   * sigma.
   * @param r result with default accuracy on output
   * @param z postion of evaluation
   * @param w postion of evaluation
   * @param k power
   */
  public void sigma(Complex r, Complex z, Complex w, int k) {
    sigma.eval(r, z, w, k, acc);
  }

  /**
   * Retunrs period matrix with prescribed accuracy.
   */
  final public ComplexMatrix getPeriodMatrix(double accuracy) {
    ComplexMatrix B = new ComplexMatrix(numGenerators);
    getPeriodMatrix(B, accuracy);
    return B;
  }

  /**
   * Retunrs period matrix with default accuracy.
   */
  final public ComplexMatrix getPeriodMatrix() {
    ComplexMatrix B = new ComplexMatrix(numGenerators);
    getPeriodMatrix(B, acc);
    return B;
  }

  /**
   * Computes period matrix with default accuracy.
   * @param B period matrix on output
   */
  public void getPeriodMatrix(ComplexMatrix B) {
    periodMatrix.eval(B, acc);
  }

  /**
   * Computes period matrix with prescribed accuracy.
   * @param B period matrix on output
   * @param accuracy of computed period matrix
   */
  public void getPeriodMatrix(ComplexMatrix B, double accuracy) {
    periodMatrix.eval(B, accuracy);
  }

  public void abelMapDifferential(ComplexVector v, Complex z) {
    v.newSize(numGenerators);

    for (int i = 0; i < numGenerators; i++) {
      abelianDifferentialOf1stKind(tmp, z, i);
      v.set(i, tmp);
    }
  }

  public void abelMap(Complex[] v, Complex z) {
    if (v.length == numGenerators) {
      for (int i = 0; i < numGenerators; i++) {
        abelianIntegralOf1stKind(v[i], z, i);
      }
    }
  }

  public void abelMap(ComplexVector v, Complex z) {
    v.newSize(numGenerators);

    for (int i = 0; i < numGenerators; i++) {
      abelianIntegralOf1stKind(tmp, z, i);
      v.set(i, tmp);
    }
  }

  /**
   * Retunrs k( sigma(j,n) ).
   */
  final double k(int j, int n, Complex P) {

    final double k = center[n][j == 0 ? 1 : 0].dist(P) - radius[n];

    return k > 0 ? k : 0;
  }

  final double k(SchottkyGroupElement sigma, Complex P) {
    return dist(sigma, P);
  }

  /**
   * Retunrs K( sigma(j,n) ).
   */
  final double K(int j, int n, Complex P) {

    return center[n][j == 0 ? 1 : 0].dist(P) + radius[n];
  }


  final double K(SchottkyGroupElement sigma, Complex P) {

     if (sigma.left == SchottkyGroupElement.IDENTITY) {
       return dist(P);
     }

     final double targetRadius = targetRadius(sigma, tmp);

     return tmp.dist(P) + targetRadius;
  }

  /** returns kappa for wordLength l. */
  public final double kappa( int l ) {
    return abelianDifferential.kappa( l );
  }



  final private double evalTheta(SchottkyGroupElement element,
                                 double theta, int wordLength) {

    if (element.wordLength > wordLength + 1)
      return theta;

    if (element.updateID != updateID) {
      updateElement(element);
    }

    if (element.wordLength == wordLength + 1) {
      theta = Math.max( theta( element ), theta );
    }

    if (element.child == null) {
      createLeftChilds(element);
    }

    final SchottkyGroupElement[] child = element.child;

    final int numOfChilds = child.length;

    for (int i = 0; i < numOfChilds; i++) {
      theta = evalTheta(child[i], theta, wordLength);
    }
    return theta;
  }

  /**
   * Computes theta for given word lengths.
   */
  public final double theta(int wordLength) {

    if (wordLength == 0 )
      return 1;

    return evalTheta( id, 0, wordLength);
  }

  /**
   * returns theta1( i, n )
   */
  final double theta1(int i, int n) {
    double max = radius[n] / k(i, n, center[n][i]);

    for (int l = 0; l < numGenerators; l++) {
      if (l != n) {
        final double t = Math.max(radius[l] / k(i, n, center[l][0]),
                                  radius[l] / k(i, n, center[l][1]));

        if (t > max) {
          max = t;
        }
      }
    }

    return max;
  }

  /**
   * retunrs maximum of theta( sigma ) with word length 1.
   */
  final public double theta1() {
    double max = 0;
    for (int i = 0; i < 2; i++) {
      for (int n = 0; n < numGenerators; n++) {
        final double theta = theta1(i, n);
        if (theta > max) {
          max = theta;
        }
      }
    }
    return max;
  }

  /**
   * retunrs maximum of theta( sigma ) with word length 1.
   */
  final public double theta2() {
    double max = 0;

    for (int i = 0; i < 2; i++) {
      for (int n = 0; n < numGenerators; n++) {

        final Complex[] innerCenter = this.innerCenter[n][i];
        final double[] innerRadius = this.innerRadius[n][i];

        for (int j = 0; j < 2; j++) {
          for (int m = 0; m < numGenerators; m++) {

            if (n != m || i != j) {
              for (int k = 0; k < innerRadius.length; k++) {
                double value
                    = radius[m] /
                    (innerCenter[k].dist(center[m][j]) - innerRadius[k]);
                if (value > max) {
                  max = value;
                }
              }
            }
          }
        }
      }
    }
    return max;
  }


  final double theta(SchottkyGroupElement tau ) {
    // throws NullPointerException for tau = id
    return radius[ tau.left ] / k( tau.parent, center[tau.left][tau.leftIsInvert] );
  }



  final double r(double q_) {
    final double q = (2 * numGenerators - 1) * q_;

    return q_ / (1 + q_) / (1 - q);
  }

  final double rMinus(double r, double q_) {
    return r - q_ / (1 - q_ * q_);
  }

  final double rMinus(double q_) {
    return r(q_) - q_ / (1 - q_ * q_);
  }

  final double rPlus(double r, double q_) {
    return r + 1 / (1 - q_ * q_);
  }

  final double rPlus(double q_) {
    return r(q_) + 1 / (1 - q_ * q_);
  }

  final Moebius leftParent = new Moebius();

  final Moebius nextLeftSubword(SchottkyGroupElement sigma) {

    leftParent.assignTimes(sigma,
                           sigma.rightIsInvert == 0 ?
                           generatorInv[sigma.right] : generator[sigma.right]);
    return leftParent;
  }

  final double targetRadius(SchottkyGroupElement sigma,
                            Complex targetCenter) {

    final double radiusOfFirst = radius[sigma.right];

    Complex centerOfFirst = center[sigma.right][sigma.rightIsInvert ==
        0 ? 1 : 0];

    if (sigma.wordLength == 1) {
      targetCenter.assign(centerOfFirst);
      return radiusOfFirst;
    }

    return nextLeftSubword(sigma).getRadiusOfMappedCircle(centerOfFirst,
        radiusOfFirst, targetCenter);
  }

  /**
   * returns dist of P to image of the fundamental domain F
   * under transformation sigma.
   * To do this really efficiently it would be best to have
   * the next left parent (word without the most right parent),
   * but our tree grows to the left, and only the right parent
   * (word without the most left parent) is given.
   * Thus we have to generate the right perent numerically.
   * @return dist( sigma(F), P )
   */
   final double dist(SchottkyGroupElement sigma, Complex P) {

     if (sigma.left == SchottkyGroupElement.IDENTITY) {
       return dist(P);
     }

     final double targetRadius = targetRadius(sigma, tmp);

     final double dist = tmp.dist(P) - targetRadius;

     if (dist > -1E-12) {
       return dist < 0 ? 0 : dist;
     }

     // P is in the image of the first circle.
     // Check if P is in one of the images of the schottky circles.
     // This part of the routine is not needed in the algorithm
     // for selecting the elements therefor we print a warning.

     System.out.println("Warning: unexpected distance measure");

     for (int j = 0; j < 2; j++) {
       for (int i = 0; i < numGenerators; i++) {

         final double theRadius = sigma.getRadiusOfMappedCircle(center[i][j],
             radius[i],
             tmp);

         final double aDist = theRadius - tmp.dist(P);

         if (aDist > 0) {
           return aDist;
         }
       }
     }

     return 0;
   }

   final Complex2By2 tmpMatrix   = new Complex2By2();
   final Moebius inverseOfParent = new Moebius();

   final double targetRadiusForInverseTransformation
       (SchottkyGroupElement sigma, Complex targetCenter) {

     final double radiusOfFirst = radius[sigma.left];

     Complex centerOfFirst = center[sigma.left][sigma.leftIsInvert ==
         0 ? 0 : 1];

     if (sigma.wordLength == 1) {
       targetCenter.assign(centerOfFirst);
       return radiusOfFirst;
     }

     tmpMatrix.assignAdjugate(sigma.parent);
     inverseOfParent.assign( tmpMatrix );

     return inverseOfParent.getRadiusOfMappedCircle(centerOfFirst,
         radiusOfFirst, targetCenter);
   }

  /** Returns distance of P to fundamental domain F.
   * Thus the result is only not zero if P lays in the interior
   * of one of the 2N circles which forms the boundary of F. */
  final double dist(Complex P) {

    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < numGenerators; i++) {

        double aDist = radius[i] - center[i][j].dist(P);

        if (aDist >= 0) {
          return aDist;
        }
      }
    }

    return 0;
  }

  final double k2(int j, int m, Complex z) {

    double result = Double.MAX_VALUE;

    final int J = j == 0 ? 1 : 0;

    Complex[] someInnerCenter = innerCenter[m][J];
    double[] someInnerRadius = innerRadius[m][J];

    for (int k = 0; k < 2 * numGenerators - 1; k++) {

      double tmpD = someInnerCenter[k].dist(z) - someInnerRadius[k];

      if (tmpD < result) {
        result = tmpD;
      }
    }

    return result;
  }


  final double d2(Complex z) {

    double result = Double.MAX_VALUE;

    double minDist = Double.MAX_VALUE;

    int minN = 0, minI = 0;

    for (int n = 0; n < numGenerators; n++) {
      for (int i = 0; i < 2; i++) {

        double thisDist = distToCenter[n][i] = center[n][i].dist(z);

        if (thisDist < minDist) {
          minDist = thisDist;
          minN = n;
          minI = i;
        }
      }
    }

    // compute dist first for nearest Circle

    Complex[] someInnerCenter = innerCenter[minN][minI];

    double[] someInnerRadius = innerRadius[minN][minI];

    for (int k = 0; k < 2 * numGenerators - 1; k++) {

      double tmpD = someInnerCenter[k].dist(z) - someInnerRadius[k];

      if (tmpD < result) {
        result = tmpD;
      }
    }

    // compute dist for all other

    for (int n = 0; n < numGenerators; n++) {
      for (int i = 0; i < 2; i++) {

        if (n != minN && i != minI && distToCenter[n][i] + radius[n] < result) {

          someInnerCenter = innerCenter[n][i];
          someInnerRadius = innerRadius[n][i];

          for (int k = 0; k < 2 * numGenerators - 1; k++) {

            double tmpD = someInnerCenter[k].dist(z) - someInnerRadius[k];

            if (tmpD < result) {
              result = tmpD;
            }
          }
        }
      }
    }

    return result;
  }

  final public void gamma(final Complex r, final double acc) {
    abelianDifferential.gamma(r, acc);
  }

  final public void gamma(final Complex r) {
    abelianDifferential.gamma(r, acc);
  }

  final public Complex gamma(final double acc) {
    final Complex r = new Complex();
    abelianDifferential.gamma(r, acc);
    return r;
  }

  final public Complex gamma() {
    final Complex r = new Complex();
    abelianDifferential.gamma(r, acc);
    return r;
  }

  final public void chi(final Complex r, final double acc) {
    abelianDifferential.chi(r, acc);
  }

  final public void chi(final Complex r) {
    abelianDifferential.chi(r, acc);
  }

  final public Complex chi(final double acc) {
    final Complex r = new Complex();
    abelianDifferential.chi(r, acc);
    return r;
  }

  final public Complex chi() {
    final Complex r = new Complex();
    abelianDifferential.chi(r, acc);
    return r;
  }
}