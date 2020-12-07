package nRalgorithm;

import java.util.ArrayList;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

public class ModifiedNR {

	public static void main(String[] args) {
		ArrayList<ArrayList<Admittance>> admittance = new ArrayList<>();
		int N=6;
		ArrayList<Bus> buses=null;
		double[] delta=null;
		double[][] theta=null;
		double [] nq = new double[N];
		double w=0;
		
		double j13[] = new double[N - 1];
		int j = 1;
		for (int k = 2; k < N; k++) {
			double temp0 = 0;
			for (int n = 1; n < N; n++) {
				double temp1 = 0;
				double admittanceX = admittance.get(k).get(n).X;
				double admittanceR = admittance.get(k).get(n).R;
				temp1 = (FastMath.pow(admittanceX, 2) / w)
						/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3 / 2));
				double temp2 = 0;
				temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
				temp0 = temp0
						+ temp1 * buses.get(n).voltage.getValue() * FastMath.cos(delta[k] - delta[n] - theta[k][n])
						+ temp2 * FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2))
								* buses.get(n).voltage.getValue() * FastMath.sin(delta[k] - delta[n] - theta[k][n]);

			}

			j13[j++] = temp0;

		}

		double j14[] = new double[N - 1];
		j = 1;
		for (int k = 2; k < N; k++) {
			double admittanceX = admittance.get(k).get(1).X;
			double admittanceR = admittance.get(k).get(1).R;
			double normAdmittance = FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2));
			j14[j++] = buses.get(k).voltage.getValue() * normAdmittance
					* FastMath.cos(delta[k] - delta[1] - theta[k][1]);
		}

		double j23[] = new double[N - 1];
		j = 1;
		for (int k = 2; k < N; k++) {
			double temp0 = 0;
			for (int n = 1; n < N; n++) {
				double temp1 = 0;
				double admittanceX = admittance.get(k).get(n).X;
				double admittanceR = admittance.get(k).get(n).R;
				temp1 = (FastMath.pow(admittanceX, 2) / w)
						/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3 / 2));
				double temp2 = 0;
				temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
				temp0 = temp0
						+ temp1 * buses.get(n).voltage.getValue() * FastMath.sin(delta[k] - delta[n] - theta[k][n])
						- temp2 * FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2))
								* buses.get(n).voltage.getValue() * FastMath.cos(delta[k] - delta[n] - theta[k][n]);

			}

			j23[j++] = temp0 * buses.get(k).voltage.getValue();

		}

		double j24[] = new double[N - 1];
		j = 1;
		for (int k = 2; k < N; k++) {
			double admittanceX = admittance.get(k).get(1).X;
			double admittanceR = admittance.get(k).get(1).R;
			double normAdmittance = FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2));
			j14[j++] = buses.get(k).voltage.getValue() * normAdmittance
					* FastMath.sin(delta[k] - delta[1] - theta[k][1]);
		}

		double j31[] = new double[N - 1];
		double j32[] = new double[N - 1];
		for (int k = 1; k < N; k++) {
			j31[k] = 0;
			j32[k] = 0;
			
		}
		double j34 = 0;
		
		double j33 = 0 ;
			
		for (Bus b : buses) {
			if (b.type == 4) 
				j33 += -1 / b.mp;
		}
		
		

		double j42[] = new double[N - 1];
		for (int k = 1; k < N; k++) {
			if(buses.get(k).type==4) {
				j42[k] = -1/nq[k];
			}else {
				j42[k] = 0;
			}
			
		}
		
		
		double j43 = 0;
		
		double j44 = 0;
		if (buses.get(1).type == 4)
			j44 = -1 / buses.get(1).nq;
		else
			j44 = 0;

			

		
		

		

		

	}

}
