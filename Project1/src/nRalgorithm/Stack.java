package nRalgorithm;

import java.util.ArrayList;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

public class Stack {
	public static double[][] constructJacabian(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, ArrayList<Bus> buses, double w,
			Complex[][] cAdmittances, ArrayList<ArrayList<Integer>> indexes) {

		int nofP = pq.get(0).size();
		int nofQ = pq.get(1).size();
		int nofD = deltaVoltageOrders.get(0).size();
		int nofV = deltaVoltageOrders.get(1).size();
		
		double[][] jacobian = new double[nofP + nofQ + 2][nofD + nofV + 2];

		for (int i = 0; i < nofP; i++) {
			for (int j = 0; j < nofD; j++) {
				jacobian[i][j] = pq.get(0).get(i)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(0).get(j)));
			}
		}

		for (int i = 0; i < nofP; i++) {
			for (int j = nofD; j < nofD + nofV; j++) {
				jacobian[i][j] = pq.get(0).get(i)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(1).get((j - nofD))));
			}
		}

		for (int i = nofP; i < nofP + nofQ; i++) {
			for (int j = 0; j < nofD; j++) {
				jacobian[i][j] = pq.get(1).get(i - nofP)
						.getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(0).get(j)));
			}
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			for (int j = nofD; j < nofD + nofV; j++) {
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

		for (int k : pS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;
				double admittanceX = cAdmittances[k][n].getImaginary();
				double admittanceR = cAdmittances[k][n].getReal();
				temp1 = (FastMath.pow(admittanceX, 2) / w)
						/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
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
			jacobian[i][nofD + nofV] = j13[i];
		}

		double j23[] = new double[nofQ];
		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(1));
		qS.addAll(indexes.get(2));
		j = 0;
		for (int k : qS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;
				double admittanceX = cAdmittances[k][n].getImaginary();
				double admittanceR = cAdmittances[k][n].getReal();
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

		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV] = j23[i - nofP];
		}
		double j33 = 0;
		for (Bus b : buses) {
			if (b.type == 2)
				j33 += -1 / b.mp;
		}
		jacobian[nofP + nofQ][nofD + nofV] = j33;

		double j43 = 0;
		jacobian[nofP + nofQ + 1][nofP + nofQ] = j43;

		double j14[] = new double[buses.size() - 1];
		j = 0;
		for (int k = 1; k < buses.size()-1; k++) {
			double magnitude = Admittance.getMag(cAdmittances[k][0]);
			j14[j++] = buses.get(k).voltage.getValue() * magnitude * FastMath.cos(
					buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - buses.get(k).theta.getEntry(k, 0));
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV + 1] = j14[i];
		}

		double j24[] = new double[nofQ];
		j = 0;

		for (int k = 1; k < nofQ; k++) {
			double magnitude = Admittance.getMag(cAdmittances[k][0]);

			j24[j] = buses.get(k).voltage.getValue() * magnitude * FastMath.sin(
					buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - buses.get(k).theta.getEntry(k, 0));
			j++;
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV + 1] = j24[i - nofP];
		}

		double j34 = 0;
		if(buses.get(0).type==2)
			jacobian[nofP + nofQ][nofD + nofV + 1] = j34;
		else
			jacobian[nofP + nofQ][nofD + nofV + 1] = -1 / (2*buses.get(0).nq);;
		

		double j44 = 0;
		if (buses.get(0).type == 2)
			j44 = -1 / (2*buses.get(0).nq);
		else
			j44 = 0;

		jacobian[nofP + nofQ + 1][nofD + nofV + 1] = j44;

		double j31 = 0;

		for (int i = 0; i < buses.size() - 1; i++) {
			jacobian[nofP + nofQ][i] = j31;

		}

		double j41 = 0;

		for (int i = 0; i < buses.size() - 1; i++) {
			jacobian[nofP + nofQ + 1][i] = j41;

		}
		double j32[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses.get(k).type == 2) {
				j32[j++] = -1 / (2*buses.get(k).nq);
			} else {
				j32[j++] = 0;
			}

		}
		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ][i] = j32[i-nofD];

		}
		
		
		double j42[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses.get(k).type == 2) {
				j42[j++] = -1 / (2*buses.get(k).nq);
			} else {
				j42[j++] = 0;
			}

		}

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ + 1][i] = j42[i - nofD];

		}

		return jacobian;

	}
	public static double[] calculateMismatchMatrix(ArrayList<Bus> Buses, double w, double w0, double v0,
			ArrayList<ArrayList<DerivativeStructure>> pq, ArrayList<Double> PQLossLoad) {

		double pSys = 0;
		double qSys = 0;
		double pTot = PQLossLoad.get(0) + PQLossLoad.get(2);
		double qTot = PQLossLoad.get(1) + PQLossLoad.get(3);

		for (Bus b : Buses) {
			if (b.type == 2 && b.index!=0) {
				pSys += b.p;
				qSys += b.q;
			}
		}

		int nofP = pq.get(0).size();
		int nofQ = pq.get(1).size();

		double[] pMismatch = new double[nofP];
		double[] qMismatch = new double[nofQ];
		for (int i = 0; i < pq.get(0).size(); i++)
			pMismatch[i] = pq.get(0).get(i).getValue();

		for (int i = 0; i < pq.get(1).size(); i++)
			qMismatch[i] = pq.get(1).get(i).getValue();

		double[] mismatches = new double[nofP + nofQ + 2];
		System.arraycopy(pMismatch, 0, mismatches, 0, pMismatch.length);
		System.arraycopy(qMismatch, 0, mismatches, pMismatch.length, qMismatch.length);
		mismatches[(pMismatch.length + qMismatch.length)] = pTot - pSys;
		mismatches[(pMismatch.length + qMismatch.length) + 1] = qTot- qSys;

//		System.out.printf("Ptotal:\t%f \nQtotal:\t %f\n", pTot, qTot);
//		System.out.printf("Psys:\t%f \nQsys:\t %f\n", pSys, qSys);
//		System.out.printf("Ptot-Psys:\t%f\n", pTot - pSys);
//		System.out.printf("Qtot-Qsys:\t%f\n", qTot - qSys);
//		System.out.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		// for(int i=0; i<mismatches.length;i++) {
		// System.out.print(mismatches[i]+"\t");
		// }
//		System.out.println();
		return mismatches;
	}
	
	
	
}
