package nRalgorithm;


import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;



public class Sketch3 {

	public static void main(String[] args) {
		int size = 5;
     NormalDistribution n = new NormalDistribution(0.1,0.01);

		
		double a[][] = new double[size][size];
		double b[][] = new double[size][size];

		fillArray(a, b);
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[1].length; j++) {
				System.out.printf("%.2f\t", a[i][j]);
			}
			System.out.println();
		}
		System.out.println();
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[1].length; j++) {
				System.out.printf("%.2f\t", b[i][j]);
			}
			System.out.println();
		}
		
		double [] arr = new double [10];
		
		
		for(int i=0 ; i<arr.length;i++) {
			arr[i]= -n.sample();
			System.out.println(arr[i]);
		}


	}

	public static void fillArray(double[][] a, double b[][]) {
		Random random = new Random();
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[1].length; j++) {
				a[i][j] = 0;
				b[i][j] = 0;
			}
		}
	
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[1].length; j++) {
				if (i < j) {
						
						int chosen = random.nextInt(a[1].length);
						if (chosen > i && findNonzeros(a[i]) < 3 ) {
							a[i][chosen] = random.nextDouble() * 0.5;
							b[i][chosen] = random.nextDouble() * 0.5;
						
					}
				} else if (i > j) {
					a[i][j] = a[j][i];
					b[i][j] = b[j][i];

				} else if (i == j) {
					b[i][j] = 0;
					a[i][j] = 0;

				}
			}
		}
		boolean flag=false;
		for(int i=0;i<a.length;i++){
			if(findNonzeros(a[i]) < 2)
				flag =true;
		}
		
		if(flag)
			fillArray(a,b);
	}

	
	public static void fillArray2(double[][] a, double b[][]) {
		Random random = new Random();
		
		for (int i = 0; i < a.length; i++) {
			int chosen = random.nextInt(1234) % a.length;
			int chosen2 = random.nextInt(1234) % a.length;
			a[i][chosen] = random.nextDouble() * 0.5;
			a[i][chosen2] = random.nextDouble() * 0.5;
		}
		
		
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[1].length; j++) {
				if (i > j) {
					a[i][j] = a[j][i];
					b[i][j] = b[j][i];

				} else {
					b[i][j] = 0;
					a[i][j] = 0;

				}
			}
		}
		boolean flag=false;
		for(int i=0;i<a.length;i++){
			if(findNonzeros(a[i]) == 0)
				flag =true;
		}
		
		if(flag)
			fillArray2(a,b);
	}
	public static int findNonzeros(double[] a) {
		int count = 0;
		for (int i = 0; i < a.length; i++)
			if(a[i]!=0)
				count++;
		return count;
	}

}
