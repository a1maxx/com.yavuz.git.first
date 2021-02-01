package nRalgorithm;



import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

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
						cAdmittances[i][j] = temp.pow(-1).multiply(-1);
					}
			}

		}
		for (int i = 0; i < radmittances.length; i++)  {
			  
			for (int j = 0; j < radmittances[1].length; j++) {
				if (i==j) {
					Complex temp0 = new Complex(0, 0);
					for (int k = 0; k < radmittances.length; k++) {
						for (int l = 0; l < radmittances[1].length; l++) {
							if (k!=l && k==i) {
								temp0 = temp0
										.add((new Complex(radmittances[k][l], xadmittances[k][l] * w)).pow(-1));
							}

						}
					}
					cAdmittances[i][j] = temp0;
				}
			}

			
		}

//		for (int i = 0; i < radmittances.length; i++) {
//			for (int j = 0; j < radmittances[1].length; j++) {
//				System.out.printf(" %.2f\t + j%.2f\t",cAdmittances[i][j].getReal(),cAdmittances[i][j].getImaginary());
//			}
//			System.out.println();
//		}
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
	
	public static double getMag(Complex x) {
		return FastMath.sqrt(Math.pow(x.getReal(),2)+Math.pow(x.getImaginary(),2));
	}
	
	public static double getAng(Complex x) {
		return FastMath.atan2(x.getImaginary(),x.getReal());
	}

	public static RealMatrix createMadmittance(Complex[][] cAdmittances) {
		double [][] mAdmittance = new double[cAdmittances.length][cAdmittances.length];
		for(int i=0;i<cAdmittances.length;i++) {
			for(int j=0;j<cAdmittances[0].length;j++) {
				mAdmittance[i][j]=Admittance.getMag(cAdmittances[i][j]);
				
			}
		}
		
		return new Array2DRowRealMatrix(mAdmittance);
	}
	
	public static RealMatrix createTadmittance(Complex[][] cAdmittances) {
		double [][] tAdmittance = new double[cAdmittances.length][cAdmittances.length];
		for(int i=0;i<cAdmittances.length;i++) {
			for(int j=0;j<cAdmittances[0].length;j++) {
				tAdmittance[i][j]=Admittance.getAng(cAdmittances[i][j]);
				
			}
		}
		
		return new Array2DRowRealMatrix(tAdmittance);
	}

	public static Complex[][] constructComplexAdmittanceMatrix2(double[][] radmittances, double[][] xadmittances,
			double w) {
		Complex[][] cAdmittances = new Complex[radmittances.length][radmittances[1].length];
		
	for (int i = 0; i < radmittances.length; i++) {
			for (int j = 0; j < radmittances[1].length; j++) {
				if (i != j && (radmittances[i][j] !=0 || xadmittances [i][j]!=0)) {
					Complex temp = new Complex(radmittances[i][j], xadmittances[i][j] * w);
					cAdmittances[i][j] = temp.pow(-1).multiply(-1);
				}else if(i!=j){
					cAdmittances[i][j] = new Complex(0,0);
				}
			}

		}
		for (int i = 0; i < radmittances.length; i++) {
			for (int j = 0; j < radmittances[1].length; j++) {
				if (i == j ) {
					Complex temp0 = new Complex(0, 0);
					for (int k = 0; k < radmittances.length; k++) {
						for (int l = 0; l < radmittances[1].length; l++) {
							if ( k!=l && k == i && (radmittances[k][l] != 0 || xadmittances [k][l]!=0)) {
								temp0 = temp0.add((new Complex(radmittances[k][l], xadmittances[k][l] * w)).pow(-1));
							}

						}
					}
					
					cAdmittances[i][j] = radmittances[i][j] * xadmittances[i][j] != 0
							? temp0.add(new Complex(radmittances[i][j], xadmittances[i][j] * w).pow(-1))
							: temp0;
				}
			}

		}

		return cAdmittances;
	}
	
	
	public void fillImpedance() {
		
		
	}
	
	
}

