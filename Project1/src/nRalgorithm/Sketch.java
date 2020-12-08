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
		
		double [] b = {1,2} ;
		double [] k = {3,4} ; 
		double [][] g = new double[2][b.length];
		System.arraycopy(b, 0, g[0], 0, b.length);
		System.arraycopy(k, 0, g[1], 0, k.length);
		RealMatrix c = new Array2DRowRealMatrix(g);
		System.out.println(c);
		
		
//		System.out.println(Arrays.asList(a).indexOf(3));
//		GenericNR.update(unknowns);
		Complex f = new Complex(1,5);
		
//		System.out.println(c.conjugate().multiply(c));
//		System.out.println(c.getReal());
//		System.out.println(c.getImaginary());
//		Complex c2= new Complex(1,2);
//		System.out.println(c2.multiply(c));
//		System.out.println(c.subtract(new Complex(0,1)));
//		
//		
		
		
	}

	

}
