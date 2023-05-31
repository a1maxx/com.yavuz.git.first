package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

public class ModifiedNR {

	public static ArrayList<ArrayList<Integer>> identifyNet(ArrayList<Bus> buses) {

		ArrayList<ArrayList<Integer>> indexes = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> PVindex = new ArrayList<Integer>();
		ArrayList<Integer> PQindex = new ArrayList<Integer>();
		ArrayList<Integer> DRindex = new ArrayList<Integer>();
		indexes.add(PVindex);
		indexes.add(PQindex);
		indexes.add(DRindex);
		for (Bus b : buses) {
			if (b.index != 0) {
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

	public static void setActiveReactiveGen(ArrayList<Bus> buses, double wi, double w0, double v0) {
		for (Bus b : buses) {
			if (b.type == 2) {
//				b.p = Math.max(Math.min(0.5*( ((1 / b.mp) * (w0 - wi)) + ((1 / b.nq) * (v0 - b.voltage.getValue())) ),b.qMax),0);		
//				b.q = Math.max(Math.min(0.5*( ((1 / b.nq) * (v0 - b.voltage.getValue())) - ((1 / b.mp) * (w0 - wi)) ),b.qMax),0);
				b.p = Math.min(0.5 * (((1 / b.mp) * (w0 - wi)) + ((1 / b.nq) * (v0 - b.voltage.getValue()))), b.pMax);
				b.q = Math.min(0.5 * (((1 / b.nq) * (v0 - b.voltage.getValue())) - ((1 / b.mp) * (w0 - wi))), b.qMax);
			}
		}

	}

	public static void setActiveReactiveGen2(ArrayList<Bus> buses, double wi, double w0, double v0) {

		for (Bus b : buses) {
			if (b.type == 2) {
				b.p = Math.min((1 / b.mp) * (w0 - wi), 10e2);
				b.q = Math.min((1 / b.nq) * (v0 - b.voltage.getValue()), 10e2);
			}
		}

	}

	public static ArrayList<Double> calculatePQLossLoad(Complex[][] cAdmittances, ArrayList<Bus> buses, double w) {

		ArrayList<Double> PQLossLoad = new ArrayList<>();

		double PLoad = 0.0;
		double QLoad = 0.0;
		double v0 = 1.0;
		double w0 = 1.0;

		double alpha = 1.0;
		double beta = 1.0;
		double kpf = 1.0;
		double kqf = -1.0;

		for (int k = 0; k < buses.size(); k++) {
			if (buses.get(k).type == 1) {
				buses.get(k).p = buses.get(k).nominal_p * Math.pow(buses.get(k).voltage.getValue() / v0, alpha)
						* (1 + kpf * (w - w0));
				PLoad += buses.get(k).p;
				buses.get(k).q = buses.get(k).nominal_q * Math.pow(buses.get(k).voltage.getValue() / v0, beta)
						* (1 + kqf * (w - w0));
				QLoad += buses.get(k).q;
			}
		}
		double PLoss = 0.0;
		double QLoss = 0.0;
		for (int k = 0; k < buses.size(); k++) {
			for (int n = 0; n < buses.size(); n++) {
			
				Complex temp0 = buses.get(k).cVolt.conjugate().multiply(buses.get(n).cVolt);
				Complex temp1 = temp0.add(buses.get(n).cVolt.conjugate().multiply(buses.get(k).cVolt));
				Complex temp2 = cAdmittances[k][n].multiply(temp1);

				PLoss = PLoss + temp2.getReal();
				QLoss = QLoss + temp2.getImaginary();
				
			}
		}

		PQLossLoad.add(0.5 * PLoss);
		PQLossLoad.add(-0.5 * QLoss);

		PQLossLoad.add(PLoad);
		PQLossLoad.add(QLoad);

		return PQLossLoad;
	}

	public static double[] calculateMismatchMatrix2(ArrayList<Bus> Buses, ArrayList<ArrayList<DerivativeStructure>> pq,
			ArrayList<Double> PQLossLoad, ArrayList<ArrayList<Integer>> indexes) {
		ArrayList<Integer> pS = new ArrayList<Integer>();
		pS.addAll(indexes.get(2));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(0));
		Collections.sort(pS);

		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(2));
		qS.addAll(indexes.get(1));
		Collections.sort(qS);

		double pSys = 0;
		double qSys = 0;
		double pTot = PQLossLoad.get(2) - PQLossLoad.get(0);
		double qTot = PQLossLoad.get(3) - PQLossLoad.get(1);

		for (Bus b : Buses) {
			if (b.type == 2 && b.index != 0) {
				pSys += b.p;
				qSys += b.q;
			} else if (b.type == 0) {
				pSys += b.p;
				qSys += b.q;
			}

		}

		int nofP = pq.get(0).size();
		int nofQ = pq.get(1).size();

		double[] pMismatch = new double[nofP];
		double[] qMismatch = new double[nofQ];
		int j = 0;
		for (int i : pS) {
			pMismatch[j] = pq.get(0).get(j).getValue() - Buses.get(i).p;
			j++;
		}
		j = 0;
		for (int i : qS) {
			qMismatch[j] = pq.get(1).get(j).getValue() - Buses.get(i).q;
			j++;
		}
		double[] mismatches = new double[nofP + nofQ + 2];
		for (int i = 0; i < nofP; i++) {
			mismatches[i] = pMismatch[i];
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			mismatches[i] = qMismatch[i - nofP];
		}

		mismatches[nofP + nofQ] = pTot + pSys;
		mismatches[nofP + nofQ + 1] = qTot + qSys;

		return mismatches;
	}

	public static double[] calculateMismatchMatrix3(ArrayList<Bus> Buses, ArrayList<ArrayList<DerivativeStructure>> pq,
			ArrayList<Double> PQLossLoad, ArrayList<ArrayList<Integer>> indexes) {
		ArrayList<Integer> pS = new ArrayList<Integer>();
		pS.addAll(indexes.get(2));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(0));

		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(2));
		qS.addAll(indexes.get(1));

		double pSys = 0;
		double qSys = 0;
		double pTot = PQLossLoad.get(2) - PQLossLoad.get(0);
		double qTot = PQLossLoad.get(3) - PQLossLoad.get(1);

		for (Bus b : Buses) {
			if (b.type == 2 && b.index != 0) {
				pSys += b.p;
				qSys += b.q;
			} else if (b.type == 0) {
				pSys += b.p;
				qSys += b.q;
			}

		}

		int nofP = pq.get(0).size();
		int nofQ = pq.get(1).size();

		double[] pMismatch = new double[nofP];
		double[] qMismatch = new double[nofQ];
		int j = 0;
		for (int i : pS) {
			pMismatch[j] = Buses.get(i).p - pq.get(0).get(j).getValue();
			j++;
		}
		j = 0;
		for (int i : qS) {
			qMismatch[j] = Buses.get(i).q - pq.get(1).get(j).getValue();
			j++;
		}

		double[] mismatches = new double[nofP + nofQ + 2];
		for (int i = 0; i < nofP; i++) {
			mismatches[i] = pMismatch[i];
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			mismatches[i] = qMismatch[i - nofP];
		}

		mismatches[nofP + nofQ] = pTot - pSys;
		mismatches[nofP + nofQ + 1] = qTot - qSys;

		return mismatches;
	}

	public static double convergenceCheck(RealMatrix R1, RealMatrix R2) {
		double absum = 0;

		for (int i = 0; i < R1.getRowDimension(); i++) {
			absum += Math.abs(R1.getEntry(i, 0) - R2.getEntry(i, 0));
		}
		return absum / R1.getRowDimension();

	}

	public static double[][] constructJacabian2(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, ArrayList<Bus> buses, double w, Complex[][] cAdmittances,
			ArrayList<ArrayList<Integer>> indexes, RealMatrix admittance, RealMatrix theta, double[][] radmittances,
			double[][] xadmittances) {

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
		pS.addAll(indexes.get(2));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(0));

		for (int k : pS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;
				double admittanceX = -cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = -cAdmittances[k][n].pow(-1).getReal();
				;
				temp1 = -(FastMath.pow(admittanceX, 2) / w)
						/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
				double temp2 = 0;
				temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
				temp0 = temp0
						+ temp1 * buses.get(n).voltage.getValue()
								* FastMath.cos(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
										- theta.getEntry(k, n))
						+ temp2 * admittance.getEntry(k, n) * buses.get(n).voltage.getValue() * FastMath.sin(
								buses.get(k).delta.getValue() - buses.get(n).delta.getValue() - theta.getEntry(k, n));

			}

			j13[j] = temp0 * buses.get(k).voltage.getValue();
			j++;

		}

		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV] = j13[i];
		}

		double j23[] = new double[nofQ];

		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(2));
		qS.addAll(indexes.get(1));

		j = 0;
		for (int k : qS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;
				double admittanceX = -cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = -cAdmittances[k][n].pow(-1).getReal();

				temp1 = -(FastMath.pow(admittanceX, 2) / w)
						/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
				double temp2 = 0;
				temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
				temp0 = temp0
						+ temp1 * buses.get(n).voltage.getValue()
								* FastMath.sin(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
										- theta.getEntry(k, n))
						- temp2 * admittance.getEntry(k, n) * buses.get(n).voltage.getValue() * FastMath.cos(
								buses.get(k).delta.getValue() - buses.get(n).delta.getValue() - theta.getEntry(k, n));

			}

			j23[j] = temp0 * buses.get(k).voltage.getValue();
			j++;

		}

		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV] = j23[i - nofP];
		}
		double j33 = 0;
		for (Bus b : buses) {
			if (b.type == 2)
				j33 += -1 / (2 * b.mp);
		}

		jacobian[nofP + nofQ][nofD + nofV] = j33;

		double j43 = 0;
		for (Bus b : buses) {
			if (b.type == 2) {
				j43 += 1 / (2 * b.mp);
			}
		}
		jacobian[nofP + nofQ + 1][nofD + nofV] = j43;

		double j14[] = new double[nofP];
		j = 0;
		for (int k : pS) {
			j14[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.cos(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV + 1] = j14[i];
		}

		double j24[] = new double[nofQ];
		j = 0;

		for (int k : qS) {

			j24[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.sin(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV + 1] = j24[i - nofP];
		}

		double j34 = 0;
		if (buses.get(0).type == 2)
			jacobian[nofP + nofQ][nofD + nofV + 1] = -1 / (2 * buses.get(0).nq);
		else
			jacobian[nofP + nofQ][nofD + nofV + 1] = j34;

		double j44 = 0;
		if (buses.get(0).type == 2)
			j44 = -1 / (2 * buses.get(0).nq);
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
				j32[j] = -1 / (2 * buses.get(k).nq);
				j++;
			} else {
				j32[j] = 0;
				j++;
			}

		}
		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ][i] = j32[i - nofD];

		}

		double j42[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses.get(k).type == 2) {
				j42[j] = -1 / (2 * buses.get(k).nq);
				j++;
			} else {
				j42[j] = 0;
				j++;
			}

		}

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ + 1][i] = j42[i - nofD];

		}

		return jacobian;

	}

	public static double[][] constructJacobian3(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, ArrayList<Bus> buses, double w, Complex[][] cAdmittances,
			ArrayList<ArrayList<Integer>> indexes, RealMatrix admittance, RealMatrix theta, double[][] radmittances,
			double[][] xadmittances) {

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

		ArrayList<Integer> pS = new ArrayList<Integer>();
		pS.addAll(indexes.get(2));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(0));
		Collections.sort(pS);

		double j13[] = new double[nofP];
		int j = 0;
		for (int k : pS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;

				double admittanceX = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();

				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 1.5));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses.get(n).voltage.getValue()
									* FastMath.cos(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n))
							+ temp2 * admittance.getEntry(k, n) * buses.get(n).voltage.getValue()
									* FastMath.sin(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}
			}

			j13[j] = temp0 * buses.get(k).voltage.getValue();
			j++;

		}

		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV] = j13[i];
		}

		double j23[] = new double[nofQ];

		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(2));
		qS.addAll(indexes.get(1));
		Collections.sort(qS);

		j = 0;
		for (int k : qS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;

				double admittanceX = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();

				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses.get(n).voltage.getValue()
									* FastMath.sin(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n))
							- temp2 * admittance.getEntry(k, n) * buses.get(n).voltage.getValue()
									* FastMath.cos(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}

			}

			j23[j] = temp0 * buses.get(k).voltage.getValue();
			j++;

		}

		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV] = j23[i - nofP];
		}
		double j33 = 0;
		for (Bus b : buses) {
			if (b.type == 2)
				j33 += -1 / (2 * b.mp);
		}

		jacobian[nofP + nofQ][nofD + nofV] = j33;

		double j43 = 0;
		for (Bus b : buses) {
			if (b.type == 2) {
				j43 += 1 / (2 * b.mp);
			}
		}
		jacobian[nofP + nofQ + 1][nofD + nofV] = j43;

		double j14[] = new double[nofP];
		j = 0;
		for (int k : pS) {
			j14[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.cos(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV + 1] = j14[i];
		}

		double j24[] = new double[nofQ];
		j = 0;

		for (int k : qS) {

			j24[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.sin(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV + 1] = j24[i - nofP];
		}

		double j34 = 0;
		if (buses.get(0).type == 2)
			jacobian[nofP + nofQ][nofD + nofV + 1] = -1 / (2 * buses.get(0).nq);
		else
			jacobian[nofP + nofQ][nofD + nofV + 1] = j34;

		double j44 = 0;
		if (buses.get(0).type == 2)
			j44 = -1 / (2 * buses.get(0).nq);
		else
			j44 = 0;

		jacobian[nofP + nofQ + 1][nofD + nofV + 1] = j44;

		double j31 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ][i] = j31;

		}

		double j41 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ + 1][i] = j41;

		}
		double j32[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses.get(k).type == 2) {
				j32[j] = -1 / (2 * buses.get(k).nq);
				j++;
			} else {
				j32[j] = 0;
				j++;
			}

		}
		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ][i] = j32[i - nofD];

		}

		double j42[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses.get(k).type == 2) {
				j42[j] = -1 / (2 * buses.get(k).nq);
				j++;
			} else {
				j42[j] = 0;
				j++;
			}

		}

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ + 1][i] = j42[i - nofD];

		}

		return jacobian;

	}

	public static double[][] constructJacobian4(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, ArrayList<Bus> buses, double w, Complex[][] cAdmittances,
			ArrayList<ArrayList<Integer>> indexes, RealMatrix admittance, RealMatrix theta, double[][] radmittances,
			double[][] xadmittances) {

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

		ArrayList<Integer> pS = new ArrayList<Integer>();
		pS.addAll(indexes.get(2));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(0));
		Collections.sort(pS);

		double j13[] = new double[nofP];
		int j = 0;
		for (int k : pS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;
				double admittanceX = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();

				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses.get(n).voltage.getValue()
									* FastMath.cos(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n))
							+ temp2 * admittance.getEntry(k, n) * buses.get(n).voltage.getValue()
									* FastMath.sin(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}
			}

			j13[j] = temp0 * buses.get(k).voltage.getValue();
			j++;

		}

		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV] = j13[i];
		}

		double j23[] = new double[nofQ];

		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(2));
		qS.addAll(indexes.get(1));
		Collections.sort(qS);

		j = 0;
		for (int k : qS) {
			double temp0 = 0;
			for (int n = 0; n < buses.size(); n++) {
				double temp1 = 0;

				double admittanceX = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();

				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses.get(n).voltage.getValue()
									* FastMath.sin(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n))
							- temp2 * admittance.getEntry(k, n) * buses.get(n).voltage.getValue()
									* FastMath.cos(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}

			}

			j23[j] = temp0 * buses.get(k).voltage.getValue();
			j++;

		}

		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV] = j23[i - nofP];
		}
		double j33 = 0;
		for (Bus b : buses) {
			if (b.type == 2)
				j33 += -1 / (b.mp);
		}

		jacobian[nofP + nofQ][nofD + nofV] = j33;

		double j43 = 0;

		jacobian[nofP + nofQ + 1][nofD + nofV] = j43;

		double j14[] = new double[nofP];
		j = 0;
		for (int k : pS) {
			j14[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.cos(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV + 1] = j14[i];
		}

		double j24[] = new double[nofQ];
		j = 0;

		for (int k : qS) {

			j24[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.sin(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV + 1] = j24[i - nofP];
		}

		double j34 = 0;
		jacobian[nofP + nofQ][nofD + nofV + 1] = j34;

		double j44 = 0;
		if (buses.get(0).type == 2)
			j44 = -1 / (buses.get(0).nq);
		else
			j44 = 0;

		jacobian[nofP + nofQ + 1][nofD + nofV + 1] = j44;

		double j31 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ][i] = j31;

		}

		double j41 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ + 1][i] = j41;

		}

		double j32[] = new double[nofV];

		for (int k = 0; k < qS.size(); k++)
			j32[k] = 0;

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ][i] = j32[i - nofD];

		}

		double j42[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses.get(k).type == 2) {
				j42[j] = -1 / (buses.get(k).nq);
				j++;
			} else {
				j42[j] = 0;
				j++;
			}

		}

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ + 1][i] = j42[i - nofD];

		}

		return jacobian;

	}

	public static double[][] artificialJacob(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders, ArrayList<Bus> buses,
			double w, ArrayList<ArrayList<Integer>> indexes, double[][] radmittances, double[][] xadmittances,
			RealMatrix DELTA, RealMatrix X0) {

		Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, w);
		RealMatrix admittance = Admittance.createMadmittance(cAdmittances);
		RealMatrix theta = Admittance.createTadmittance(cAdmittances);

		double alpha = 0.0;
		double beta = 0.0;
		double kpf = 1.0;
		double kqf = -1.0;

		double w0 = 1.00;
		double v0 = 1.01;

		int params = buses.size() * 2;
		int order = 2;

		RealMatrix X1 = X0.add(DELTA.scalarMultiply(0.5));

		@SuppressWarnings("unchecked")
		ArrayList<Bus> buses2 = (ArrayList<Bus>) buses.clone();

		ModifiedNR.updateUnknowns(X1, buses2, deltaVoltageOrders, params, order);
		w = X1.getEntry(X1.getRowDimension() - 2, 0);
		buses2.get(0).voltage = new DerivativeStructure(params, 2, 1, X1.getEntry(X1.getRowDimension() - 1, 0));

		ModifiedNR.setActiveReactiveGen(buses2, w, w0, v0);

		for (int k = 0; k < buses2.size(); k++) {
			if (buses2.get(k).type == 1) {
				buses2.get(k).p = buses2.get(k).nominal_p * Math.pow(buses2.get(k).voltage.getValue() / v0, alpha)
						* (1 + kpf * (w - w0));

				buses2.get(k).q = buses2.get(k).nominal_q * Math.pow(buses2.get(k).voltage.getValue() / v0, beta)
						* (1 + kqf * (w - w0));

			}
		}

		ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses2, admittance, theta);

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
		ArrayList<Integer> pS = new ArrayList<Integer>();
		pS.addAll(indexes.get(2));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(0));
		Collections.sort(pS);

		int j = 0;
		for (int k : pS) {
			double temp0 = 0;
			for (int n = 0; n < buses2.size(); n++) {
				double temp1 = 0;
				double admittanceX = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();

//						&& admittanceR !=0 && admittanceX != 0
				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses2.get(n).voltage.getValue()
									* FastMath.cos(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n))
							+ temp2 * admittance.getEntry(k, n) * buses2.get(n).voltage.getValue()
									* FastMath.sin(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}
			}

			j13[j] = temp0 * buses2.get(k).voltage.getValue();
			j++;

		}

		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV] = j13[i];
		}

		double j23[] = new double[nofQ];

		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(2));
		qS.addAll(indexes.get(1));
		Collections.sort(qS);

		j = 0;
		for (int k : qS) {
			double temp0 = 0;
			for (int n = 0; n < buses2.size(); n++) {
				double temp1 = 0;

				double admittanceX = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();
//						
				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses2.get(n).voltage.getValue()
									* FastMath.sin(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n))
							- temp2 * admittance.getEntry(k, n) * buses2.get(n).voltage.getValue()
									* FastMath.cos(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}

			}

			j23[j] = temp0 * buses2.get(k).voltage.getValue();
			j++;

		}

		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV] = j23[i - nofP];
		}
		double j33 = 0;
		for (Bus b : buses2) {
			if (b.type == 2)
				j33 += -1 / (2 * b.mp);
		}

		jacobian[nofP + nofQ][nofD + nofV] = j33;

		double j43 = 0;
		for (Bus b : buses2) {
			if (b.type == 2) {
				j43 += 1 / (2 * b.mp);
			}
		}
		jacobian[nofP + nofQ + 1][nofD + nofV] = j43;

		double j14[] = new double[nofP];
		j = 0;
		for (int k : pS) {
			j14[j] = buses2.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.cos(buses2.get(k).delta.getValue() - buses2.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV + 1] = j14[i];
		}

		double j24[] = new double[nofQ];
		j = 0;

		for (int k : qS) {

			j24[j] = buses2.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.sin(buses2.get(k).delta.getValue() - buses2.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV + 1] = j24[i - nofP];
		}

		double j34 = 0;
		if (buses2.get(0).type == 2)
			jacobian[nofP + nofQ][nofD + nofV + 1] = -1 / (2 * buses2.get(0).nq);
		else
			jacobian[nofP + nofQ][nofD + nofV + 1] = j34;

		double j44 = 0;
		if (buses2.get(0).type == 2)
			j44 = -1 / (2 * buses2.get(0).nq);
		else
			j44 = 0;

		jacobian[nofP + nofQ + 1][nofD + nofV + 1] = j44;

		double j31 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ][i] = j31;

		}

		double j41 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ + 1][i] = j41;

		}
		double j32[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses2.get(k).type == 2) {
				j32[j] = -1 / (2 * buses2.get(k).nq);
				j++;
			} else {
				j32[j] = 0;
				j++;
			}

		}
		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ][i] = j32[i - nofD];

		}

		double j42[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses2.get(k).type == 2) {
				j42[j] = -1 / (2 * buses2.get(k).nq);
				j++;
			} else {
				j42[j] = 0;
				j++;
			}

		}

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ + 1][i] = j42[i - nofD];

		}

		return jacobian;

	}

	public static double[][] artificialJacob2(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders, ArrayList<Bus> buses,
			double w, ArrayList<ArrayList<Integer>> indexes, double[][] radmittances, double[][] xadmittances,
			RealMatrix DELTA, RealMatrix X0) {

		Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, w);
		RealMatrix admittance = Admittance.createMadmittance(cAdmittances);
		RealMatrix theta = Admittance.createTadmittance(cAdmittances);

		double alpha = 0.0;
		double beta = 0.0;
		double kpf = 1.0;
		double kqf = -1.0;

		double w0 = 1.00;
		double v0 = 1.01;

		int params = buses.size() * 2;
		int order = 2;

		RealMatrix X1 = X0.add(DELTA.scalarMultiply(0.5));

		@SuppressWarnings("unchecked")
		ArrayList<Bus> buses2 = (ArrayList<Bus>) buses.clone();

		ModifiedNR.updateUnknowns(X1, buses2, deltaVoltageOrders, params, order);
		w = X1.getEntry(X1.getRowDimension() - 2, 0);
		buses2.get(0).voltage = new DerivativeStructure(params, 2, 1, X1.getEntry(X1.getRowDimension() - 1, 0));

		ModifiedNR.setActiveReactiveGen2(buses2, w, w0, v0);

		for (int k = 0; k < buses2.size(); k++) {
			if (buses2.get(k).type == 1) {
				buses2.get(k).p = buses2.get(k).nominal_p * Math.pow(buses2.get(k).voltage.getValue() / v0, alpha)
						* (1 + kpf * (w - w0));

				buses2.get(k).q = buses2.get(k).nominal_q * Math.pow(buses2.get(k).voltage.getValue() / v0, beta)
						* (1 + kqf * (w - w0));

			}
		}

		ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses2, admittance, theta);

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
		ArrayList<Integer> pS = new ArrayList<Integer>();
		pS.addAll(indexes.get(2));
		pS.addAll(indexes.get(1));
		pS.addAll(indexes.get(0));
		Collections.sort(pS);

		int j = 0;
		for (int k : pS) {
			double temp0 = 0;
			for (int n = 0; n < buses2.size(); n++) {
				double temp1 = 0;
				double admittanceX = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k != n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();

				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses2.get(n).voltage.getValue()
									* FastMath.cos(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n))
							+ temp2 * admittance.getEntry(k, n) * buses2.get(n).voltage.getValue()
									* FastMath.sin(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}
			}

			j13[j] = temp0 * buses2.get(k).voltage.getValue();
			j++;

		}

		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV] = j13[i];
		}

		double j23[] = new double[nofQ];

		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(2));
		qS.addAll(indexes.get(1));
		Collections.sort(qS);

		j = 0;
		for (int k : qS) {
			double temp0 = 0;
			for (int n = 0; n < buses2.size(); n++) {
				double temp1 = 0;

				double admittanceX = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getImaginary()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getImaginary();
				double admittanceR = k == n
						? cAdmittances[k][n].pow(-1).isNaN() ? 0 : -cAdmittances[k][n].pow(-1).getReal()
						: cAdmittances[k][n].pow(-1).isNaN() ? 0 : cAdmittances[k][n].pow(-1).getReal();
//						
				if (admittance.getEntry(k, n) != 0) {
					temp1 = -(FastMath.pow(admittanceX, 2) / w)
							/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
					double temp2 = 0;
					temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
					temp0 = temp0
							+ temp1 * buses2.get(n).voltage.getValue()
									* FastMath.sin(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n))
							- temp2 * admittance.getEntry(k, n) * buses2.get(n).voltage.getValue()
									* FastMath.cos(buses2.get(k).delta.getValue() - buses2.get(n).delta.getValue()
											- theta.getEntry(k, n));
				}

			}

			j23[j] = temp0 * buses2.get(k).voltage.getValue();
			j++;

		}

		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV] = j23[i - nofP];
		}

		double j33 = 0;
		for (Bus b : buses) {
			if (b.type == 2)
				j33 += -1 / (b.mp);
		}

		jacobian[nofP + nofQ][nofD + nofV] = j33;

		double j43 = 0;

		jacobian[nofP + nofQ + 1][nofD + nofV] = j43;

		double j14[] = new double[nofP];
		j = 0;
		for (int k : pS) {
			j14[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.cos(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV + 1] = j14[i];
		}

		double j24[] = new double[nofQ];
		j = 0;

		for (int k : qS) {

			j24[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath
					.sin(buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV + 1] = j24[i - nofP];
		}

		double j34 = 0;
		jacobian[nofP + nofQ][nofD + nofV + 1] = j34;

		double j44 = 0;
		if (buses.get(0).type == 2)
			j44 = -1 / (buses.get(0).nq);
		else
			j44 = 0;

		jacobian[nofP + nofQ + 1][nofD + nofV + 1] = j44;

		double j31 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ][i] = j31;

		}

		double j41 = 0;

		for (int i = 0; i < nofD; i++) {
			jacobian[nofP + nofQ + 1][i] = j41;

		}

		double j32[] = new double[nofV];

		for (int k = 0; k < qS.size(); k++)
			j32[k] = 0;

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ][i] = j32[i - nofD];

		}

		double j42[] = new double[nofV];
		j = 0;
		for (int k : qS) {
			if (buses.get(k).type == 2) {
				j42[j] = -1 / (buses.get(k).nq);
				j++;
			} else {
				j42[j] = 0;
				j++;
			}

		}

		for (int i = nofD; i < nofD + nofV; i++) {
			jacobian[nofP + nofQ + 1][i] = j42[i - nofD];

		}

		return jacobian;
	}

	public static ArrayList<ArrayList<Integer[]>> createOrders2(ArrayList<Bus> buses) {

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = new ArrayList<>();
		ArrayList<Integer[]> deltaOrders = new ArrayList<Integer[]>();
		ArrayList<Integer[]> voltageOrders = new ArrayList<Integer[]>();

		for (Bus b : buses) {
			if (b.index != 0) {
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
					if (buses.get(k).type != 0) {
						equationP = equationP.add(buses.get(k).voltage.multiply(buses.get(k).admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
										.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).cos()));

						equationQ = equationQ.add(buses.get(k).voltage.multiply(-buses.get(k).admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
										.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).sin()));
					} else {
						equationP = equationP.add(buses.get(k).voltage.multiply(buses.get(k).admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
										.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).cos()));

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

	public static ArrayList<ArrayList<DerivativeStructure>> createEquations4(ArrayList<Bus> buses,
			RealMatrix admittance, RealMatrix theta) {
		ArrayList<ArrayList<DerivativeStructure>> pq = new ArrayList<ArrayList<DerivativeStructure>>();
		ArrayList<DerivativeStructure> p = new ArrayList<DerivativeStructure>();
		ArrayList<DerivativeStructure> q = new ArrayList<DerivativeStructure>();
		for (int k = 1; k < buses.size(); k++) {
			DerivativeStructure equationP = new DerivativeStructure(buses.size() * 2, 2);
			DerivativeStructure equationQ = new DerivativeStructure(buses.size() * 2, 2);
			for (int i = 0; i < buses.size(); i++) {
				if (admittance.getEntry(k, i) != 0) {
					if (buses.get(k).type != 0) {
						equationP = equationP.add(buses.get(k).voltage.multiply(admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply(
										(buses.get(k).delta.negate().add(theta.getEntry(k, i)).add(buses.get(i).delta))
												.cos()));

						equationQ = equationQ.add(buses.get(k).voltage.multiply(-admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply(
										(buses.get(k).delta.negate().add(theta.getEntry(k, i)).add(buses.get(i).delta))
												.sin()));
					} else {
						equationP = equationP.add(buses.get(k).voltage.multiply(admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply(
										(buses.get(k).delta.negate().add(theta.getEntry(k, i)).add(buses.get(i).delta))
												.cos()));

					}
				}

			}
			p.add(equationP);
			if (buses.get(k).type != 0) {
				q.add(equationQ);
			}
		}
		pq.add(p);
		pq.add(q);

		return pq;
	}

	public static ArrayList<ArrayList<DerivativeStructure>> createEquations5(ArrayList<Bus> buses,
			RealMatrix admittance, RealMatrix theta) {
		ArrayList<ArrayList<DerivativeStructure>> pq = new ArrayList<ArrayList<DerivativeStructure>>();
		ArrayList<DerivativeStructure> p = new ArrayList<DerivativeStructure>();
		ArrayList<DerivativeStructure> q = new ArrayList<DerivativeStructure>();
		for (int k = 1; k < buses.size(); k++) {
			DerivativeStructure equationP = new DerivativeStructure(buses.size() * 2, 2);
			DerivativeStructure equationQ = new DerivativeStructure(buses.size() * 2, 2);
			for (int i = 0; i < buses.size(); i++) {
				if (admittance.getEntry(k, i) != 0) {
					if (buses.get(k).type != 0) {
						equationP = equationP.add(buses.get(i).voltage.multiply(admittance.getEntry(k, i)).multiply(
								(buses.get(k).delta.add(buses.get(i).delta.negate()).add(-theta.getEntry(k, i)))
										.cos()));

						equationQ = equationQ.add(buses.get(i).voltage.multiply(admittance.getEntry(k, i)).multiply(
								(buses.get(k).delta.add(buses.get(i).delta.negate()).add(-theta.getEntry(k, i)))
										.sin()));

					} else {
						equationP = equationP.add(buses.get(i).voltage.multiply(admittance.getEntry(k, i)).multiply(
								(buses.get(k).delta.add(-theta.getEntry(k, i)).add(buses.get(i).delta.negate()))
										.cos()));
					}
				}
			}
			p.add(equationP.multiply(buses.get(k).voltage));
			if (buses.get(k).type != 0) {
				q.add(equationQ.multiply(buses.get(k).voltage));
			}
		}
		pq.add(p);
		pq.add(q);

		return pq;
	}

	public static RealMatrix setUnknownValues(ArrayList<Bus> buses, ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			double w, double v1) {
		int nofD = deltaVoltageOrders.get(0).size();
		int nofV = deltaVoltageOrders.get(1).size();
		double[] deltas = new double[nofD];
		double[] voltages = new double[nofV];
		double[] unknowns = new double[nofD + nofV + 2];
		int i = 0;
		int j = 0;
		for (Bus b : buses) {
			if (b.index != 0) {
				if (b.type == 0) {
					deltas[i] = b.delta.getValue();
					i++;
				} else if (b.type == 1) {
					deltas[i] = b.delta.getValue();
					i++;
					voltages[j] = b.voltage.getValue();
					j++;
				} else if (b.type == 2) {
					deltas[i] = b.delta.getValue();
					i++;
					voltages[j] = b.voltage.getValue();
					j++;
				} else {
					System.out.println("ERROR!");
				}
			}

		}

		System.arraycopy(deltas, 0, unknowns, 0, nofD);
		System.arraycopy(voltages, 0, unknowns, nofD, nofV);
		unknowns[nofD + nofV] = w;
		unknowns[nofD + nofV + 1] = v1;
		RealMatrix X0 = new Array2DRowRealMatrix(unknowns);

		return X0;
	}

	public static RealMatrix setUnknownValues2(ArrayList<Bus> buses, ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			double w, double v1) {
		int nofD = deltaVoltageOrders.get(0).size();
		int nofV = deltaVoltageOrders.get(1).size();
		double[] deltas = new double[nofD];
		double[] voltages = new double[nofV];
		double[] unknowns = new double[nofD + nofV + 2];
		int i = 0;
		int j = 0;
		for (Bus b : buses) {
			if (b.index != 0) {
				if (b.type == 0) {
					deltas[i] = b.delta.getValue();
					i++;
				} else if (b.type == 1 || b.type == 2) {
					deltas[i] = b.delta.getValue();
					i++;
					voltages[j] = b.voltage.getValue();
					j++;
				} else {
					System.out.println("ERROR!");
				}
			}

		}

		System.arraycopy(deltas, 0, unknowns, 0, nofD);
		System.arraycopy(voltages, 0, unknowns, nofD, nofV);
		unknowns[nofD + nofV] = w;
		unknowns[nofD + nofV + 1] = v1;
		RealMatrix X0 = new Array2DRowRealMatrix(unknowns);

		return X0;
	}

	public static void updateUnknowns(RealMatrix X1, ArrayList<Bus> buses,
			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders, int params, int order) {
		int nofD = deltaVoltageOrders.get(0).size();
		int j = 0;
		int f = nofD;
		for (int i = 0; i < deltaVoltageOrders.size(); i++) {
			for (Integer[] o : deltaVoltageOrders.get(i)) {
				int k = Arrays.asList(o).indexOf(1);
				if (k % 2 == 0) {
					buses.get(k / 2).delta = new DerivativeStructure(params, order, k, X1.getEntry(j, 0));
					j++;
				} else if (k % 2 == 1) {
					buses.get(k / 2).voltage = new DerivativeStructure(params, order, k, X1.getEntry(f, 0));
					f++;
				}

			}

		}

		for (Bus b : buses) {
			b.cVolt = Admittance.polarToComplex(b.voltage.getValue(), b.delta.getValue());
		}

	}

	public static double sumMatrix(RealMatrix x) {
		double sum = 0;
		for (int i = 0; i < x.getRowDimension(); i++) {
			for (int j = 0; j < x.getColumnDimension(); j++) {
				sum += Math.abs(x.getEntry(i, j));
			}
		}
		return sum;

	}
	
	public static double evaluate(ArrayList<Bus> buses,ArrayList<Double> pqLossLoad, double wi,double [] position,RealMatrix fx0) {

		
		double fitness = 50*buses.size()*(10*Math.abs(pqLossLoad.get(0)) + Math.abs(pqLossLoad.get(1)));
		
		for (Bus b : buses) {
			fitness = fitness + Math.abs(1 - b.voltage.getValue());

		}

		for(int i =0 ; i<position.length;i++) {
			if(position[i]<0)
				fitness -= 100*position[i];
//			else if(position[i]>1)
//				fitness += 100*position[i];
										
		}
		
		if (ModifiedNR.sumMatrix(fx0) > 1E-4)
			fitness += ModifiedNR.sumMatrix(fx0);

		fitness += (Math.abs(1 - wi) * 10) * buses.size();

		
		return fitness;
	

		
	}
	
	public static double evaluate2(ArrayList<Bus> buses,ArrayList<Double> pqLossLoad, double wi,double [] position,RealMatrix fx0) {

		
		double fitness = 200*(10*Math.abs(pqLossLoad.get(0)) + Math.abs(pqLossLoad.get(1)));
		
		for (Bus b : buses) {
			if(b.voltage.getValue()<0.95)
				fitness += 100*Math.abs(0.95-b.voltage.getValue());
			else if(b.voltage.getValue()>1.05)
				fitness += 100* Math.abs(b.voltage.getValue()-1.05);

		}
		
		if(wi<0.995)
			fitness += 1000*Math.abs(0.995-wi);
		else if(wi>1.005)
			fitness += 1000*Math.abs(wi-1.005);
		
		for(int i =0 ; i<position.length;i++) {
			if(position[i]<0)
				fitness -= 100*position[i];
			else if(position[i]>1)
				fitness += 100*position[i];
										
		}
		
		if (ModifiedNR.sumMatrix(fx0) > 1E-4)
			fitness += ModifiedNR.sumMatrix(fx0);

	
		
		return fitness;
	

		
	}

	public static double evaluate3(ArrayList<Bus> buses,ArrayList<Double> pqLossLoad, double wi,double [] position,RealMatrix fx0) {

		ArrayList<Double> vals = new ArrayList<Double>();

//		double fitness = 200*(10*Math.abs(pqLossLoad.get(0)) + Math.abs(pqLossLoad.get(1)));
		
		double fitness; 
		
		for (Bus b : buses) {
			if(b.voltage.getValue()<0.95)
				vals.add(Math.abs(0.95-b.voltage.getValue()));
			else if(b.voltage.getValue()>1.05)
				vals.add(Math.abs(b.voltage.getValue()-1.05));

		}
		
		if(wi<0.995)
			vals.add(Math.abs(0.995-wi));
		else if(wi>1.005)
			vals.add(Math.abs(wi-1.005));
		
//		for(int i =0 ; i<position.length;i++) {
//			if(position[i]<0)
//				vals.add(-position[i]);
//			else if(position[i]>1)
//				vals.add(position[i]-1);									
//		}
		
		
		if(vals.isEmpty()) 
			fitness = Math.atan(10*Math.abs(pqLossLoad.get(0)) + Math.abs(pqLossLoad.get(1))) - Math.PI/2;
		else
			fitness = Collections.max(vals);
		
		
		if (ModifiedNR.sumMatrix(fx0) > 1E-3)
			fitness += ModifiedNR.sumMatrix(fx0)*1e2;

	
		
		return fitness;
	
		
	}
	
	public static boolean checkAdequacy(ArrayList<Bus> buses, double SGMAX) {

		double demA = 0;
		double demR = 0;

		for (Bus b : buses) {
			if (b.type == 1) {
				demA += b.nominal_p;
				demR += b.nominal_q;
			}

		}
		double demT = Math.sqrt(Math.pow(demA, 2) + Math.pow(demR, 2));
		demT += demT * 0.05;
		return demT > SGMAX;

	}

}