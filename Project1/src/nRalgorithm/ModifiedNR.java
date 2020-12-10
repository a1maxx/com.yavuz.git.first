package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

public class ModifiedNR {

	public static void main(String[] args) {
				
		double Y11, Y12, Y13, Y21, Y22, Y23, Y31, Y32, Y33;
		Y11 = -14;
		Y12 = 10;
		Y13 = 4;
		Y21 = 10;
		Y22 = 15;
		Y23 = 5;
		Y31 = 4;
		Y32 = 5;
		Y33 = 5;

		double pi = Math.PI;
		double t11, t12, t13, t21, t22, t23, t31, t32, t33;
		t11 = -pi / 2;
		t12 = pi / 2;
		t13 = pi / 2;
		t21 = pi / 2;
		t22 = -pi / 2;
		t23 = pi / 2;
		t31 = pi / 2;
		t32 = pi / 2;
		t33 = -pi / 2;

		// Known variables

		int params = 6;
		int order = 2;

		RealMatrix X0 = new Array2DRowRealMatrix();
		double[][] admittances = { { Y11, Y12, Y13}, { Y21, Y22, Y23 }, { Y31, Y32, Y33 } };
		double[][] theta = { { t11, t12, t13 }, { t21, t22, t23 }, { t31, t32, t33 } };

		ArrayList<ArrayList<Admittance>> admittanceS = new ArrayList<>();

		for (int i = 0; i < 3; i++) {
			ArrayList<Admittance> temp0 = new ArrayList<>();
			for (int j = 0; j < 3; j++) {
				Admittance temp1 = new Admittance(4, 15);
				temp0.add(temp1);
			}
			admittanceS.add(temp0);
		}
		
		ArrayList<Bus> buses = new ArrayList<Bus>();
		buses.add(new Bus(2, admittances, theta, 0, 1.0, 0.95));
		buses.add(new Bus(2, admittances, theta, 2, -0.9, -0.5));
		buses.add(new Bus(0, admittances, theta, 4, 0.6, 0.2));
		buses.get(0).voltage = new DerivativeStructure(params, order, 1, 1.0);
		buses.get(0).delta = new DerivativeStructure(params, order, 0, 0.87);
		buses.get(2).voltage = new DerivativeStructure(params, order, 5, 1.01);
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = createOrders2(buses);
		ArrayList<ArrayList<DerivativeStructure>> pq = createEquations3(buses);
		int N = 3;
		double w = 1;
		for(ArrayList<Integer[]> i : deltaVoltageOrders) {
			for(int j=0;j<i.size();j++) {
				for(int k=0;k<i.get(j).length;k++) {
					System.out.print(i.get(j)[k]);
				}
				System.out.println();
			}	
		}
		
		System.out.println(deltaVoltageOrders.get(0).size());
		System.out.println(deltaVoltageOrders.get(1).size());
		System.out.println(pq.get(0).size());
		System.out.println(pq.get(1).size());
		ArrayList<ArrayList<Integer>> indexes =identifyNet(buses);
		
		double[][] Jacob = constructJacabian(deltaVoltageOrders, pq, N, buses, w, admittanceS,indexes);
		RealMatrix JJ = new Array2DRowRealMatrix(Jacob);
		System.out.println(JJ);
		 System.out.println(MatrixUtils.inverse(JJ));
//		double[] mismatches = calculateMismatchMatrix(buses, w, 1, 1, pq);

	}

	public static ArrayList<ArrayList<Integer>> identifyNet(ArrayList<Bus> buses) {

		ArrayList<ArrayList<Integer>> indexes = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> PVindex = new ArrayList<Integer>();
		ArrayList<Integer> PQindex = new ArrayList<Integer>();
		ArrayList<Integer> DRindex = new ArrayList<Integer>();
		indexes.add(PVindex);
		indexes.add(PQindex);
		indexes.add(DRindex);
		for (Bus b : buses) {
			if (b.index!=0) {
				if (b.type == 0) {
					indexes.get(0).add(buses.indexOf(b));
				} else if (b.type == 1) {
					indexes.get(1).add(buses.indexOf(b));
				} else if (b.type == 2) {
					indexes.get(2).add(buses.indexOf(b));
				} else
					System.out.println("ERROR!");
			}
		}

		return indexes;
	}

	public void setActiveReactive(ArrayList<Bus> buses, double wi, double w0, double v0) {
		for (Bus b : buses) {
			if (b.type == 2) {
				b.p = 1 / b.mp * (w0 - wi);
				b.q = 1 / b.nq * (v0 - b.voltage.getValue());
			}
		}

	}

	public static double[] calculateMismatchMatrix(ArrayList<Bus> Buses, double w, double w0, double v0,
			ArrayList<ArrayList<DerivativeStructure>> pq) {

		double pSys = 0;
		double qSys = 0;
		double pTot = 0;
		double qTot = 0;
		for (Bus b : Buses) {
			pTot += b.p;
			qTot += b.q;
			if (b.type == 2) {
				pSys += (1 / b.mp) * w0 - w;
				pSys += (1 / b.mp) * v0 - b.voltage.getValue();
			}

		}
		double[] pMismatch = new double[Buses.size()];
		double[] qMismatch = new double[Buses.size()];
		for (int i = 0; i < pq.size(); i++) {
			pMismatch[i] = -pq.get(0).get(i).getValue();
			qMismatch[i] = pq.get(1).get(i).getValue();

		}
		double[] mismatches = new double[(Buses.size()) * 2];
		System.arraycopy(pMismatch, 0, mismatches, 0, pMismatch.length);
		System.arraycopy(qMismatch, 0, mismatches, pMismatch.length, qMismatch.length);
		mismatches[(pMismatch.length + qMismatch.length) - 2] = pTot - pSys;
		mismatches[(pMismatch.length + qMismatch.length) - 1] = qTot - qSys;

		return mismatches;
	}

	public static double[][] constructJacabian(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, int N, ArrayList<Bus> buses, double w,
			ArrayList<ArrayList<Admittance>> admittance,ArrayList<ArrayList<Integer>> indexes) {
		
		int nofP =pq.get(0).size();
		int nofQ =pq.get(1).size();
		int nofD =deltaVoltageOrders.get(0).size();
		int nofV = deltaVoltageOrders.get(1).size();
		N = nofP+ nofQ;
		double[][] jacobian = new double[nofP+nofQ+2][nofD+nofV+2];
		
		
		for (int i = 0; i < nofP; i++) {
			for (int j = 0; j < nofD; j++) {
				jacobian[i][j] = pq.get(0).get(i)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(0).get(j)));
			}
		}

		for (int i = 0; i < nofP; i++) {
			for (int j = nofD; j < nofD+nofV; j++) {
				jacobian[i][j] = pq.get(0).get(i)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(1).get((j - nofD))));
			}
		}

		for (int i = nofP; i < nofP+nofQ; i++) {
			for (int j = 0; j < nofD; j++) {
				jacobian[i][j] = pq.get(1).get(i - nofP)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(0).get(j)));
			}
		}
		for (int i = nofP; i < nofP+nofQ; i++) {
			for (int j = nofD; j < nofD+nofV; j++) {
				jacobian[i][j] = pq.get(1).get(i - nofP)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(1).get(j - nofD)));
			}
		}

		
		double j13[] = new double[nofP];
		int j = 0;
		ArrayList<Integer> pS = new ArrayList<Integer>();
		pS.addAll(indexes.get(0));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(2));
		
		for (int k :pS) {
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
						+ temp1 * buses.get(n).voltage.getValue()
								* FastMath.cos(buses.get(k).delta.getValue() - buses.get(k).delta.getValue()
										- buses.get(k).theta.getEntry(k, n))
						+ temp2 * FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2))
								* buses.get(n).voltage.getValue() * FastMath.sin(buses.get(k).delta.getValue()
										- buses.get(k).delta.getValue() - buses.get(k).theta.getEntry(k, n));

			}

			j13[j++] = temp0;

		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD+nofV] = j13[i];
		}
		double j23[] = new double[N - 1];
		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(1));
		qS.addAll(indexes.get(2));
		j = 1;
		for (int k :qS) {
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
						+ temp1 * buses.get(n).voltage.getValue()
								* FastMath.sin(buses.get(k).delta.getValue() - buses.get(k).delta.getValue()
										- buses.get(k).theta.getEntry(k, n))
						- temp2 * FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2))
								* buses.get(n).voltage.getValue() * FastMath.cos(buses.get(k).delta.getValue()
										- buses.get(k).delta.getValue() - buses.get(k).theta.getEntry(k, n));

			}

			j23[j++] = temp0 * buses.get(k).voltage.getValue();

		}

		for (int i = nofP; i < nofP+nofQ; i++) {
			jacobian[i][nofD+nofV-1] = j23[i - nofP];
		}
		double j33 = 0;
		for (Bus b : buses) {
			if (b.type == 2)
				j33 += -1 / b.mp;
		}
		jacobian[nofP+nofQ][nofD+nofV] = j33;

		double j43 = 0;
		jacobian[nofP+nofQ+1][nofP+nofQ] = j43;

		double j14[] = new double[N - 1];
		j = 0;
		for (int k = 1; k < N; k++) {
			double admittanceX = admittance.get(k).get(1).X;
			double admittanceR = admittance.get(k).get(1).R;
			double normAdmittance = FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2));
			j14[j++] = buses.get(k).voltage.getValue() * normAdmittance * FastMath.cos(
					buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - buses.get(k).theta.getEntry(k, 1));
		}

		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD+nofV+1] = j14[i];
		}
		double j24[] = new double[nofQ];
		j = 0;
		
		for (int k = 1; k < nofQ; k++) {
			double admittanceX = admittance.get(k).get(1).X;
			double admittanceR = admittance.get(k).get(1).R;
			double normAdmittance = FastMath.sqrt(Math.pow(admittanceX, 2) + Math.pow(admittanceR, 2));
			j24[j] = buses.get(k).voltage.getValue() * normAdmittance * FastMath.sin(
					buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - buses.get(k).theta.getEntry(k, 1));
			j++;
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD+nofV+1] = j14[i];
		}

		double j34 = 0;
		jacobian[nofP+nofQ][nofD+nofV+1] = j34;

		double j44 = 0;
		if (buses.get(1).type == 2)
			j44 = -1 / buses.get(1).nq;
		else
			j44 = 0;

		jacobian[nofP+nofQ+1][nofD+nofV+1] = j44;

		double j31 = 0;

		for (int i = 0; i < N - 1; i++) {
			jacobian[nofP+nofQ][i] = j31;

		}

		double j41 = 0;

		for (int i = 0; i < N - 1; i++) {
			jacobian[nofP+nofQ+1][i] = j41;

		}
		double j32 = 0;

		for (int i = nofD; i < nofD+nofV; i++) {
			jacobian[nofP+nofQ][i] = j32;

		}
		double j42[] = new double[nofV];
		j=0;
		for (int k :qS){
			if (buses.get(k).type == 2) {
				j42[j++] = -1 / buses.get(k).nq;
			} else {
				j42[j++] = 0;
			}

		}

		for (int i = nofD; i < nofD+nofV ; i++) {
			jacobian[nofP+nofQ+1][i] = j42[i-nofD];

		}

		return jacobian;

	}

	public static ArrayList<ArrayList<Integer[]>> createOrders2(ArrayList<Bus> buses) {

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = new ArrayList<>();
		ArrayList<Integer[]> deltaOrders = new ArrayList<Integer[]>();
		ArrayList<Integer[]> voltageOrders = new ArrayList<Integer[]>();

		for (Bus b : buses) {
			if (b.index!=0) {
				Integer[] temp = new Integer[buses.size() * 2];
				Arrays.fill(temp, 0);
				temp[b.index] = 1;
				deltaOrders.add(temp);
				if (b.type != 0) {
					temp = new Integer[buses.size() * 2];
					Arrays.fill(temp, 0);
					temp[b.index + 1] = 1;
					voltageOrders.add(temp);
				} 
			}
		}
		deltaVoltageOrders.add(deltaOrders);
		deltaVoltageOrders.add(voltageOrders);

		return deltaVoltageOrders;
	}

	public static ArrayList<ArrayList<DerivativeStructure>> createEquations3(ArrayList<Bus> buses) {
		ArrayList<ArrayList<DerivativeStructure>> pq = new ArrayList<ArrayList<DerivativeStructure>>();
		ArrayList<DerivativeStructure> p = new ArrayList<DerivativeStructure>();
		ArrayList<DerivativeStructure> q = new ArrayList<DerivativeStructure>();
		for (int k = 1; k < buses.size(); k++) {
			DerivativeStructure equationP = new DerivativeStructure(buses.size() * 2, 2);
			DerivativeStructure equationQ = new DerivativeStructure(buses.size() * 2, 2);
			for (int i = 0; i < buses.size(); i++) {
				if (buses.get(k).admittance.getEntry(k, i) != 0) {
					equationP = equationP.add(buses.get(k).voltage.multiply(buses.get(k).admittance.getEntry(k, i))
							.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
									.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).cos()));
					if(buses.get(k).type!=0) {
					equationQ = equationQ.add(buses.get(k).voltage.multiply(-buses.get(k).admittance.getEntry(k, i))
							.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
									.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).sin()));
					}

				}

			}

			p.add(equationP.subtract(buses.get(k).p));
			if (buses.get(k).type != 0) {
				q.add(equationQ.subtract(buses.get(k).q));
			}
		}
		
		pq.add(p);
		pq.add(q);
		return pq;
	}

}
