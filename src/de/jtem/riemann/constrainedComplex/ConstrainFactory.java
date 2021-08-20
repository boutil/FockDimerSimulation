package de.jtem.riemann.constrainedComplex;

import de.jtem.riemann.surface.SurfacePoint;

public class ConstrainFactory {

	public final static int REALSYMMETRY= 0;
	public final static int IMAGINARYSYMMETRY= 1;
	public final static int UNITCIRCLESYMMETRY= 2;
	public final static int FIX= 3;
	public final static int DEPENDEND= 4;
	public final static int REAL= 5;
	public final static int POINTSYMMETRY= 6;
	public final static int IMAGINARY= 7;

	private ConstrainFactory(){}
	
	public static void setFix(SurfacePoint p){
		new Fix(p);
	}
	
	public static void setImaginarySymmetry(SurfacePoint p1, SurfacePoint p2){
		new ImaginarySymmetry(p1,p2);
	}
	
	public static void setRealSymmetry(SurfacePoint p1, SurfacePoint p2){
		new RealSymmetry(p1,p2);
	}
	
	public static void setUnitCircleSymmetry(SurfacePoint p1, SurfacePoint p2){
		new UnitCircleSymmetry(p1,p2);
	}

	public static void setDependend(SurfacePoint p){
		new Dependent(p);
	}
	
	public static void setReal(SurfacePoint p) {
		new Real(p);
	}

	public static void setImaginary(SurfacePoint p) {
		new Imaginary(p);
	}
	
	public static void setPointSymmetry(SurfacePoint p1, SurfacePoint p2) {
		new PointSymmetry(p1,p2);
	}
	
	public static int getConstrainTypeOf(Constrain c){
		if(c instanceof Fix){
			return FIX;
		}
		if(c instanceof RealSymmetry){
			return REALSYMMETRY;
		}
		if(c instanceof ImaginarySymmetry){
			return IMAGINARYSYMMETRY;
		}
		if(c instanceof Real){
			return REAL;
		}
		if(c instanceof UnitCircleSymmetry){
			return UNITCIRCLESYMMETRY;
		}
		if(c instanceof Dependent){
			return DEPENDEND;
		}
		if(c instanceof PointSymmetry){
			return POINTSYMMETRY;
		}
		if(c instanceof Imaginary){
			return IMAGINARY;
		}
		return -1;
	}
	
	public static ConstrainedComplex[] getConstrainedComplexOf(Constrain c){
		if(c instanceof SingleConstrain){
			return new ConstrainedComplex[]{((SingleConstrain)c).A};
		}
		if(c instanceof ToupleConstrain){
			return new ConstrainedComplex[]{((ToupleConstrain)c).A,((ToupleConstrain)c).B};
		}
		return null;
	}
	
	public static boolean isFixConstrained(SurfacePoint p){
		return p.getConstrain() instanceof Fix;
	}
	
	public static boolean isImaginarySymmetryConstrained(SurfacePoint p){
		return p.getConstrain() instanceof ImaginarySymmetry;
	}
	
	public static boolean isRealSymmetryConstrained(SurfacePoint p){
		return p.getConstrain() instanceof RealSymmetry;
	}
	
	public static boolean isUnitCircleSymmetryConstrained(SurfacePoint p) {
		return p.getConstrain() instanceof UnitCircleSymmetry;
	}
	
	public static boolean isDependendConstrained(SurfacePoint p){
		return p.getConstrain() instanceof Dependent;
	}
	
	public static boolean isRealConstrained(SurfacePoint p){
		return p.getConstrain() instanceof Real;
	}
	
	public static boolean isPointSymmetryConstrained(SurfacePoint p){
		return p.getConstrain() instanceof PointSymmetry;
	}
	
	public static boolean isImaginaryConstrained(SurfacePoint p){
		return p.getConstrain() instanceof Imaginary;
	}
	
}
