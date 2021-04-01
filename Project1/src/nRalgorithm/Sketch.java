package nRalgorithm;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;


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
//		double [] d= {1,2};
//		double [] f= {3,4};
//		double [] g  = new double[5];
//		System.arraycopy(d, 0, g, 0, d.length);
//		System.arraycopy(f, 0, g, d.length, f.length);
//		for(double i : g)
//			System.out.println(i);
//		System.out.println(c.conjugate().multiply(c));
//		System.out.println(c.getReal());
//		System.out.println(c.getImaginary());
//		Complex c2= new Complex(5,2);
//		Complex c3 = new Complex(0,0.250000);
//		System.out.println(c3.pow(-1).getReal()>0.0000000001);
//		System.out.println(Math.atan2(c2.getImaginary(),c2.getReal()));
//		System.out.println(3.0/2);
//		System.out.println(0.751e-3);
//		System.out.println(2e-3*1000);
//		System.out.println(c2.multiply(c));
//		System.out.println(c.subtract(new Complex(0,1)));
//		
//		for (int i = 0; i < 3; i++) {
//			for (int j = 0; j < 3; j++) {
//				if(i!=j)
//				System.out.printf("(%d,%d)\t", i, j);
//				else System.out.print("     \t");
//			}
//			System.out.printf("\n");
//		}
//		System.out.println(Admittance.polarToComplex(1.01238, -1.47170*Math.PI/180).getReal());
		double[][] xadmittances = new double[6][6];
		double[][] theta = new double[6][6];
		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		
		File file = new File("/Users/my_mac/git/com.yavuz.git.first/Project1/test.txt");
		try {
			Scanner scan = new Scanner(file);
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					xadmittances[i][j] = scan.nextDouble();
					
				}
			}
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					radmittances[i][j]= scan.nextDouble();
					theta[i][j] = Admittance.getAng(new Complex(radmittances[i][j],xadmittances[i][j]));
				}
			}
			
//			for (int i = 0; i < xadmittances.length; i++) {
//				for (int j = 0; j < xadmittances[1].length; j++) {
//					System.out.print(radmittances[i][j]+ " " + xadmittances[i][j]+"\t");
//				}
//				System.out.println();
//			}			
			scan.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, 1);
		for (int i = 0; i < xadmittances.length; i++) {
			for (int j = 0; j < xadmittances[1].length; j++) {
				System.out.print(cAdmittances[i][j]);
			}
			System.out.println();
		}		
		
		Random random= new Random();
		
		for (int i = 0; i < xadmittances.length; i++) {
			for (int j = 0; j < xadmittances[1].length; j++) {
				xadmittances[i][j] = i == j ? 0 : (random.nextDouble() + 0.5);
			}
		}

		for (int i = 0, len = radmittances.length; i < len; i++)
			Arrays.fill(radmittances[i], (random.nextDouble()+1)*0.2);		
		
		for (int i = 0; i < xadmittances.length; i++) {
			for (int j = 0; j < xadmittances[1].length; j++) {
				theta[i][j] = Admittance.getAng(new Complex(radmittances[i][j],xadmittances[i][j]));
			}
		}
		System.out.println("---------------------------------------------------------------------------------------");
		cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, 1);
		for (int i = 0; i < xadmittances.length; i++) {
			for (int j = 0; j < xadmittances[1].length; j++) {
				System.out.print(cAdmittances[i][j]);
			}
			System.out.println();
		}		

		
	}
	

	

}
