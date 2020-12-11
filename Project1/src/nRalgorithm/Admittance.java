package nRalgorithm;

import java.util.ArrayList;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealMatrix;

public class Admittance {
	
	double X;
	double R;
	
	public Admittance(double X,double R){
		
		this.X = X;
		this.R = R;

	}
	
	public static Complex[][] constructComplexAdmittanceMatrix(double [][] radmittances,double [][] xadmittances) {
		
		Complex[][] cAdmittances = new Complex[radmittances.length][radmittances[1].length];
		for(int i=0;i<radmittances.length;i++) {
			for(int j=0;j<radmittances[1].length;j++){
				cAdmittances[i][j] = new Complex(radmittances[i][j],xadmittances[i][j]);	
			}
			
		}
		
		
		return cAdmittances;
	}

}
