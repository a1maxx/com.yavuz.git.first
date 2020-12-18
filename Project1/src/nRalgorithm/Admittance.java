package nRalgorithm;

import java.util.ArrayList;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealMatrix;

public class Admittance {
	
	double X;
	double R;
	double theta;
	double r;
	public Admittance(double X,double R){
		
		this.X = X;
		this.R = R;

	}

	public static Complex[][] constructComplexAdmittanceMatrix(double[][] radmittances, double[][] xadmittances,
			double w) {

		Complex[][] cAdmittances = new Complex[radmittances.length][radmittances[1].length];
		for (int i = 0; i < radmittances.length; i++) {
			for (int j = 0; j < radmittances[1].length; j++) {
				 	if (i!=j) {
						Complex temp = new Complex(radmittances[i][j], xadmittances[i][j] * w);
						cAdmittances[i][j] = temp.pow(-1);
					}
			}

		}
		for (int i = 0; i < radmittances.length; i++) {
			for (int j = 0; j < radmittances[1].length; j++) {
				if (i==j) {
					Complex temp0 = new Complex(0, 0);
					for (int k = 0; k < radmittances.length; k++) {
						for (int l = 0; l < radmittances[1].length; l++) {
							if (k!=l) {
								temp0 = temp0
										.add((new Complex(radmittances[k][l], xadmittances[k][l] * w)).pow(-1));
							}

						}
					}
					cAdmittances[i][j] = temp0;
				}
			}

			
		}

		for (int i = 0; i < radmittances.length; i++) {
			for (int j = 0; j < radmittances[1].length && i==j; j++) {
				System.out.print(cAdmittances[i][j].getReal());
			}
			System.out.println();
		}
		return cAdmittances;
	}

	public void getPolar(Complex C1) {
		double tmp = Math.pow(C1.getReal(), 2.0) + Math.pow(C1.getImaginary(), 2.0);
		this.r = Math.sqrt(tmp);
		this.theta = Math.atan2(C1.getImaginary(), C1.getReal());
	}
	
	public static Complex polarToComplex(double r,double theta) {
		Complex temp = new Complex(r*Math.cos(theta),r*Math.sin(theta));
		return temp;
		
	}
}

