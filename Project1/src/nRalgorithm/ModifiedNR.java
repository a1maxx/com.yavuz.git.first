package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

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
			j24[j++] = buses.get(k).voltage.getValue() * normAdmittance
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
	
	public static double[] calculateMismatchMatrix(ArrayList<Bus> Buses, double w, double w0, double v0,
			ArrayList<ArrayList<DerivativeStructure>> pq) {

		RealMatrix mismatchrMatrix = new Array2DRowRealMatrix();
		double pSys = 0;
		double qSys = 0;
		double pTot = 0;
		double qTot = 0;
		for (Bus b : Buses) {
			pTot += b.p;
			qTot += b.q;
			if (b.type == 4) {
				pSys += (1 / b.mp) * w0 - w;
				pSys += (1 / b.mp) * v0 - b.voltage.getValue();
			}

		}
		double [] pMismatch = new double[Buses.size()];
		double [] qMismatch = new double[Buses.size()];
		for(int i = 0 ; i<pq.size();i++) {
			pMismatch[i]=-pq.get(0).get(i).getValue();
			qMismatch[i]= pq.get(1).get(i).getValue();

		}
		double [] mismatches = new double[(Buses.size())*2];
		 System.arraycopy(pMismatch, 0, mismatches, 0, pMismatch.length);
		 System.arraycopy(qMismatch, 0, mismatches, pMismatch.length, qMismatch.length);
		 mismatches[pMismatch.length+qMismatch.length] = pTot-pSys;
		 mismatches[pMismatch.length+qMismatch.length+1] = qTot-qSys;
		 
		
		 return mismatches;
	}
	
	public static double[][] constructJacabian(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, int N,
			ArrayList<Bus> buses,double[] delta,
			double[][] theta,double w,ArrayList<ArrayList<Admittance>> admittance) {
		
		double[][] jacobian = new double[2 * N][2 * N];
		for (int i = 0; i < N - 1; i++) {
			for (int j = 0; j < N - 1; j++) {
				jacobian[i][j] = pq.get(0).get(i)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(0).get(j)));
			}
		}
		for (int i = N - 1; i < 2 * N - 2; i++) {
			for (int j = N - 1; j < 2 * N - 2; j++) {
				jacobian[i][j] = pq.get(0).get(i)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(1).get(j)));
			}
		}
		double j13[] = new double[N - 1];
		int j = 0;
		for (int k = 1; k < N; k++) {
			double temp0 = 0;
			for (int n = 0; n < N; n++) {
				double temp1 = 0;
				double admittanceX = admittance.get(k).get(n).X;
				double admittanceR = admittance.get(k).get(n).R;
				temp1 = (FastMath.pow(admittanceX, 2) / w)
						/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3 / 2));
				double temp2 = 0;
				temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
				temp0 = temp0
						+ temp1 * buses.get(n).voltage.getValue() * FastMath.cos(buses.get(k).delta.getValue() - buses.get(k).delta.getValue() 
								- buses.get(k).theta.getEntry(k,n))
						+ temp2 * FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2))
								* buses.get(n).voltage.getValue() * FastMath.sin(buses.get(k).delta.getValue() - buses.get(k).delta.getValue() 
										- buses.get(k).theta.getEntry(k,n));

			}

			j13[j++] = temp0;

		}
		for(int i = 0;i<N-1;i++){
			jacobian[i][2*N-2] = j13[i];	
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
		
		for(int i = N-1;i<2*N-2;i++){
			jacobian[i][2*N-1] = j23[i];	
		}
		double j33 = 0 ;
		jacobian[2*N-2][2*N-2] = j33;
		
		double j43 = 0;
		jacobian[2*N-1][2*N-2] = j43;
		
		
		double j14[] = new double[N - 1];
		j = 0;
		for (int k = 1; k < N; k++) {
			double admittanceX = admittance.get(k).get(1).X;
			double admittanceR = admittance.get(k).get(1).R;
			double normAdmittance = FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2));
			j14[j++] = buses.get(k).voltage.getValue() * normAdmittance
					* FastMath.cos(delta[k] - delta[1] - theta[k][1]);
		}
		
		for(int i = 0;i<N-1;i++){
			jacobian[i][2*N-1] = j14[i];	
		}
		double j24[] = new double[N - 1];
		j = 0;
		for (int k = 1; k < N; k++) {
			double admittanceX = admittance.get(k).get(1).X;
			double admittanceR = admittance.get(k).get(1).R;
			double normAdmittance = FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2));
			j24[j++] = buses.get(k).voltage.getValue() * normAdmittance
					* FastMath.sin(delta[k] - delta[1] - theta[k][1]);
		}
		for(int i = N-1;i<2*N-2;i++){
			jacobian[i][2*N-1] = j14[i];	
		}
		
		double j34 = 0;
		jacobian[2*N-2][2*N-1] = j34;
		
		double j44 = 0;
		if (buses.get(1).type == 4)
			j44 = -1 / buses.get(1).nq;
		else
			j44 = 0;

		jacobian[2*N-2][2*N-1] = j44;
		
		
		

				
		
		
		
		double j31[] = new double[N - 1];
		double j32[] = new double[N - 1];
		Arrays.fill(j31,0);
		Arrays.fill(j32,0);
	
			
		for (Bus b : buses) {
			if (b.type == 4) 
				j33 += -1 / b.mp;
		}
		double j42[] = new double[N - 1];
		for (int k = 1; k < N; k++) {
			if(buses.get(k).type==4) {
				j42[k] = -1/buses.get(k).nq;
			}else {
				j42[k] = 0;
			}
			
		}
	
		

		
		return jacobian;

	}
	
	static ArrayList<ArrayList<Integer[]>> createOrders2(ArrayList<Bus> buses) {

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = new ArrayList<>();
		ArrayList<Integer[]> deltaOrders = new ArrayList<Integer[]>();
		ArrayList<Integer[]> voltageOrders = new ArrayList<Integer[]>();

		for (Bus b : buses) {
				Integer[] temp = new Integer[buses.size() * 2];
				Arrays.fill(temp, 0);
				temp[b.index] = 1;
				deltaOrders.add(temp);
				temp = new Integer[buses.size() * 2];
				Arrays.fill(temp, 0);
				temp[b.index + 1] = 1;
				voltageOrders.add(temp);		
		}
		deltaVoltageOrders.add(deltaOrders);
		deltaVoltageOrders.add(voltageOrders);
		
		return deltaVoltageOrders;
	}

}
