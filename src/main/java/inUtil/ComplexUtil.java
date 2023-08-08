package inUtil;

import de.jtem.mfc.field.Complex;

public class ComplexUtil {
    public static Complex compMod(Complex x, double y){
        // x mod y behaving the same way as Math.floorMod but with doubles
        return new Complex(x.re - Math.floor(x.re/y) * y, x.im - Math.floor(x.im/y) * y);
    }
}
