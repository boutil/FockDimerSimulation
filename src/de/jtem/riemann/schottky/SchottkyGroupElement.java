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

import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;

// we define a class for matrices with left,right,norm for above iteration:
class SchottkyGroupElement
    extends Moebius
    implements Comparable {

  private static final long serialVersionUID = 1L;

  protected static final char IDENTITY = Short.MAX_VALUE;

  int updateID;

  public int left;
  public int leftIsInvert;
  public int right;
  public int rightIsInvert;

  public int wordLength;

  public SchottkyGroupElement parent;

  SchottkyGroupElement[] child;

  public double norm;

  public static long instanceCount = 0;

  Complex[] imageOfA;
  Complex[] imageOfB;

  double [] LOfInverse;

  SchottkyGroupElement() {
    instanceCount++;

    left = right = leftIsInvert = rightIsInvert = IDENTITY;

    wordLength = 0;
    parent = null;
    norm = 1.0;
  }

  public void assignInvert(SchottkyGroupElement sigma) {

    assignAdjugate(sigma);

    int tmp = sigma.left;
    left = sigma.right;
    right = tmp;

    tmp = sigma.leftIsInvert;
    leftIsInvert = sigma.rightIsInvert == 0 ? 1 : 0;
    rightIsInvert = tmp == 0 ? 1 : 0;

    parent = null;

    wordLength = sigma.wordLength;
  }

  public int compareTo(Object aTransform) {

    SchottkyGroupElement other = (SchottkyGroupElement) aTransform;

    if (other.norm < norm) {
      return 1;
    }
    else if (other.norm > norm) {
      return -1;
    }
    return 0;
  }

  void createLeftChilds(int numGenerators) {

    if (parent == null) {
      throw new RuntimeException("method does not create childs of id");
    }

    final int numberOfChilds = numGenerators * 2 - 1;

    if (child == null || child.length != numberOfChilds) {
      child = new SchottkyGroupElement[numberOfChilds];

    }
    int j = 0;

    for (char k = 0; k < left; k++) {
      child[j++] = createLeftChild(0, k);
      child[j++] = createLeftChild(1, k);
    }

    child[j++] = createLeftChild(leftIsInvert == 0 ? 0 : 1, left);

    for (char k = (char) (left + 1); k < numGenerators; k++) {
      child[j++] = createLeftChild(0, k);
      child[j++] = createLeftChild(1, k);
    }
  }

  SchottkyGroupElement createLeftChild(int i, int n) {
    if (left == n && leftIsInvert != i) {
      throw new IllegalArgumentException("does not create new word");
    }

    final SchottkyGroupElement leftChild = new SchottkyGroupElement();

    leftChild.right = right;
    leftChild.rightIsInvert = rightIsInvert;
    leftChild.parent = this;
    leftChild.left = (char) n;
    leftChild.leftIsInvert = (char) i;
    leftChild.wordLength = (char) (wordLength + 1);
    leftChild.updateID = -1;
    return leftChild;
  }

  /**
   * Computes chi = (1-2b*c) / d^4
   * @param chi on
   */
  final void getChi( final Complex chi ) {
    // Using chi temp. to compute d^4.
    chi.assign( dRe, dIm );
    chi.assignSqr();
    chi.assignSqr();

    final double d4Re = chi.re;
    final double d4Im = chi.im;

    chi.assignTimes( -2*bRe, -2*bIm, cRe, cIm );
    chi.re+=1;
    chi.assignDivide( d4Re, d4Im );
  }

  /**
   * Computes c^-2
   * @param cSqrInvert equals c^-2 on output
   */
  final void getInverseOfCSqr(final Complex inverseOfCSqr) {

    final double abs = cRe * cRe + cIm * cIm;
    final double absSqr = abs * abs;

    inverseOfCSqr.re = (cRe * cRe - cIm * cIm) / absSqr;
    inverseOfCSqr.im = -2 * cRe * cIm / absSqr;
  }

  /**
   * Computes sigma( z ) - sigma( w ), which is returned in r.
   * @param z
   * @param r = ( z-w ) / ( (cz+d) * (cw+d) )
   */
  final void diff(final Complex z, final Complex w, final Complex r) {

    final double nr = z.re - w.re;
    final double ni = z.im - w.im;

    final double dzr = cRe * z.re - cIm * z.im + dRe;
    final double dzi = cRe * z.im + cIm * z.re + dIm;

    final double dwr = cRe * w.re - cIm * w.im + dRe;
    final double dwi = cRe * w.im + cIm * w.re + dIm;

    final double dzwr = dzr * dwr - dzi * dwi;
    final double dzwi = dzr * dwi + dzi * dwr;

    final double dzw = dzwr * dzwr + dzwi * dzwi;

    r.re = (nr * dzwr + ni * dzwi) / dzw;
    r.im = (ni * dzwr - nr * dzwi) / dzw;
  }

  /**
   * Computes sigma( z ) - sigma( infty ), which is returned in r.
   * @param z
   * @param r = - 1 / ( (cz+d) * c )
   */
  final void diffWithInfty(final Complex z, final Complex r) {

    final double d1r = cRe * z.re - cIm * z.im + dRe;
    final double d1i = cRe * z.im + cIm * z.re + dIm;

    final double d2r = d1r * cRe - d1i * cIm;
    final double d2i = d1r * cIm + d1i * cRe;

    final double d2 = d2r * d2r + d2i * d2i;

    r.re = -d2r / d2;
    r.im = d2r / d2;
  }

  final Complex sz = new Complex(); // sigma( z )
  final Complex sw = new Complex(); // sigma( w );

  final Complex szsw = new Complex(); // sz * sw
  final Complex sz2 = new Complex(); // sz^2
  final Complex sw2 = new Complex(); // sw^2
  final Complex szm = new Complex(); // sz^m , m = n - 1;
  final Complex swm = new Complex(); // sw^m
  final Complex sm = new Complex(); // sz^m + sw^m
  final Complex s1 = new Complex(); // sz + sw
  final Complex d1 = new Complex(); // sz - sw
  final Complex dn = new Complex(); // sz^n - sw^n
  final Complex fm = new Complex(); // (sz-sw) sm


  Complex P = new Complex();
  Complex Q = new Complex(1);
  Complex lastQ = new Complex(1);
  
  /**
   * Computes sigma(z)^k - sigma(w)^k, which is returned in r.
   * @param z
   * @paran w
   * @param k is power
   * @param r = sigma(z)^k - sigma(w)^k on output
   */
  final void diffPow(final Complex z, final Complex w,
                     final int k, final Complex r) {

    diff(z, w, d1);
    diffPow(z, w, k, d1, r);
  }

  /**
   * Computes sigma(z)^k - sigma(w)^k, which is returned in r.
   * @param z
   * @paran w
   * @param k is power
   * @param d = element(z) - element(w), is input parameter
   * @param r = sigma(z)^k - sigma(w)^k on output
   */
  final void diffPow(final Complex z, final Complex w,
                     final int k, final Complex d, final Complex r) {

  	P.assign(0);
  	Q.assign(1);
  	
  	applyTo(z, sz);
  	applyTo(w, sw);
  	
  	szsw.assignTimes( sz, sw );
  	s1.assignPlus( sz, sw );
  	
  	for( int i=1;i<k;i++) {
  		lastQ.assign(Q);
  		P.assignTimes( szsw );
  		Q.assignTimes(s1);
  		Q.assignMinus( P );
  		P.assign( lastQ );
  	}
  	
  	r.assignTimes( d, Q );
  	/*
    if( k == 1) { 
    	r.assign(d);
    	return;
    }

    applyTo(z, sz);
    applyTo(w, sw);
    
    d1.assign( d );
    
   	if( k == 2 ) {
   		sm.assignPlus( sz, sw );
   		r.assignTimes( d1, sm );
   		return;
   		
   	} else if( k == 3 ) {
    	szm.assignSqr(sz);
    	swm.assignSqr(sw);
    	
    	szsw.assignTimes( sz, sw );
    	
    	r.assignPlus( szm, swm );
    	r.assignPlus( szsw );
    	r.assignTimes( d1 );
    	
    	return;
    	
    } else {
    	szm.assignPow( sz, k );
    	swm.assignPow( sw, k );
    	
    	r.assignMinus( szm, swm );
    	return;
    }
*/

  /*
   * 
    int n;
   
    if (k % 2 == 0) {
      n = 2; // m = 1

      applyTo(z, sz);
      applyTo(w, sw);

      dn.assignPlus(sz, sw);
      dn.assignTimes(d);

      if (k == 2) {
        r.assign(dn);
        return;
      }

      szm.assign(sz);
      swm.assign(sw);

    }
    else {
      n = 1; // m=0
      dn.assign(d);

      if (k == 1) {
        r.assign(dn);
        return;
      }

      applyTo(z, sz);
      applyTo(w, sw);

      szm.assign(1);
      swm.assign(1);
    }

    sz2.assignTimes(sz, sz);
    sw2.assignTimes(sw, sw);

    szsw.assignTimes(sz, sw);

    for (n += 2; n <= k; n += 2) {
      szm.assignTimes(sz2);
      swm.assignTimes(sw2);

      sm.assignPlus(szm, swm);
      fm.assignTimes(sm, d);

      dn.assignTimes(szsw);
      dn.assignPlus(fm);
    }

    r.assign(dn);
    
    */
  }

  /* the following functions will be removed soon */

  private static final Complex[] v = new Complex[2];
  private static final Complex[] w = new Complex[2];

  static {
    v[0] = new Complex();
    v[1] = new Complex();
  }

  static {
    w[0] = new Complex();
    w[1] = new Complex();
  }

  final public double getCircles(final Complex[] c) {
    return getCircles(c[0], c[1]);
  }

  final public double getCircles(final Complex c, final Complex c_) {
    if (isLoxodromic()) {
      double d = cRe * cRe + cIm * cIm;
      c.assign( ( -dRe * cRe - dIm * cIm) / d, (dRe * cIm - dIm * cRe) / d);
      c_.assign( (aRe * cRe + aIm * cIm) / d, (aIm * cRe - aRe * cIm) / d);
      return 1 / Math.sqrt(d);
    }
    else {
      throw new IllegalArgumentException("Transformation is not loxodromic");
    }
  }

  /**
   * Tests whether generator generate a classical Schottky group.
   * @return true iff this generate a classical Schottky group
   */
  final public boolean isClassical() {
    double d = 2 * getCircles(v);

    return Complex.distSqr(v[0], v[1]) > d * d;
  }

  /**
   * Tests wheter isometric cirlces of this and m intersect.
   */
  final boolean intersects(SchottkyGroupElement m) {
    double d = getCircles(w) + m.getCircles(v); // Attention: v is used getCircles via isLox ..
    return
        Complex.dist(v[0], w[0]) <= d ||
        Complex.dist(v[0], w[1]) <= d ||
        Complex.dist(v[1], w[0]) <= d ||
        Complex.dist(v[1], w[1]) <= d;
  }

  /**
   * Tests whether array of generators generate a classical Schottky group.
   * @param g array of generators
   * @return true iff generators generate classical Schottky group
   */
  final public static boolean isClassical(final SchottkyGroupElement[]
                                          generator) {
    final int vol = generator.length;
    for (int i = 0; i < vol; i++) {
      if (!generator[i].isClassical()) {
        return false;
      }
      for (int j = i + 1; j < vol; j++) {
        if (generator[i].intersects(generator[j])) {
          return false;
        }
      }
    }
    return true;
  }


  /**
   * Returns word generating this group element.
   * a, b, c, ... are the letters for the first, secound, thrird, ...
   * generator; Capital letters stand for their inverse.
   * @return word generating this group element.
   */

  final public String word() {
    final char [] word = new char [wordLength];

    SchottkyGroupElement e = this;

    for( int i=0; i<wordLength; i++ ) {

      word[i] = e.firstLetter();
    }

    return new String( word);
  }

  /**
   * Returns first letter of word generating this group element.
   * a, b, c, ... are the letters for the first, secound, thrird, ...
   * generator; Capital letters stand for their inverse.
   * @return first letter of word generating this group element.
   */
  final public char firstLetter() {
    return this.leftIsInvert == 1 ? (char)('A' + left) : (char)('a' + left);
  }

  final SchottkyGroupElement childWithFirstLetter( char letter ) {
    if( child == null )
      throw new RuntimeException( "no childs" );

    for( int i=0; i<child.length; i++ ) {
      if( child[i].firstLetter() == letter )
        return child[i];
    }

    throw new IllegalArgumentException( "no child has this as first letter." );
  }
}