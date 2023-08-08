package inUtil;

import java.math.BigDecimal;
import java.math.MathContext;

import de.jtem.mfc.field.Complex;

public class ComplexHighPrecision {
    
    public BigDecimal re;
    public BigDecimal im;

    public ComplexHighPrecision(Complex z) {
        re = new BigDecimal(z.re, MathContext.DECIMAL128);
        im = new BigDecimal(z.im, MathContext.DECIMAL128);
    }

    public ComplexHighPrecision(BigDecimal re, BigDecimal im) {
        this.re = re;
        this.im = im;
    }

    public ComplexHighPrecision conjugate() {
        return new ComplexHighPrecision(re, im.negate(MathContext.DECIMAL128));
    }

    public ComplexHighPrecision times(ComplexHighPrecision z) {
        return new ComplexHighPrecision(re.multiply(z.re, MathContext.DECIMAL128).subtract(im.multiply(z.im, MathContext.DECIMAL128), MathContext.DECIMAL128), im.multiply(z.re,MathContext.DECIMAL128).add(re.multiply(z.im, MathContext.DECIMAL128), MathContext.DECIMAL128));
    }

}
