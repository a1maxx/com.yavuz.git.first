package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class Sketch {

	public static void main(String[] args) {
//		ArrayList<Double> unknowns = new ArrayList<Double>();
//		if (true) {
//			Double temp = 1.0;
//			unknowns.add(temp);
//			temp = 2.0;
//			unknowns.add(temp);
//			System.out.println(unknowns);
//		}
//		RealMatrix X0 = new Array2DRowRealMatrix();
//		X0 = new Array2DRowRealMatrix(GenericNR.convertArray(unknowns));
//		System.out.println(X0);
//		System.out.println(unknowns.indexOf(2.0));
//		Integer [] a = {1,2,3};
		
//		double [] b = {1,2} ;
//		double [] k = {3,4} ; 
//		double [][] g = new double[2][b.length];
//		System.arraycopy(b, 0, g[0], 0, b.length);
//		System.arraycopy(k, 0, g[1], 0, k.length);
//		RealMatrix c = new Array2DRowRealMatrix(g);
//		System.out.println(c);
		
		
//		System.out.println(Arrays.asList(a).indexOf(3));
//		GenericNR.update(unknowns);
//		double[][] a = new double[1][2];
//		System.out.println(a.length);
//		System.out.println(a[0].length);
		double [] d= {1,2};
		double [] f= {3,4};
		double [] g  = new double[5];
		System.arraycopy(d, 0, g, 0, d.length);
		System.arraycopy(f, 0, g, d.length, f.length);
		for(double i : g)
			System.out.println(i);
//		System.out.println(c.conjugate().multiply(c));
//		System.out.println(c.getReal());
//		System.out.println(c.getImaginary());
		Complex c2= new Complex(1,2);
		System.out.println(Math.atan2(c2.getReal(), c2.getImaginary()));
		System.out.println(0.751e-3);
		System.out.println(2e-3*1000);
//		System.out.println(c2.multiply(c));
//		System.out.println(c.subtract(new Complex(0,1)));
//		
//		
		
		
	}

	

}
