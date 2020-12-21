package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

public class ModifiedNR {

	public static void main(String[] args) {

		double X11, X12, X13, X21, X22, X23, X31, X32, X33;
		X11 = 0.02;
		X12 = 0.1;
		X13 = 0.25;
		X21 = 0.1;
		X22 = 0;
		X23 = 0.2;
		X31 = 0.25;
		X32 = 0.2;
		X33 = 0.02;

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

		double[][] xadmittances = { { X11, X12, X13 }, { X21, X22, X23 }, { X31, X32, X33 } };
		double[][] theta = { { t11, t12, t13 }, { t21, t22, t23 }, { t31, t32, t33 } };

		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		for (int i = 0, len = radmittances.length; i < len; i++)
			Arrays.fill(radmittances[i], 0.04);

		ArrayList<Bus> buses = new ArrayList<Bus>();

		buses.add(new Bus(0, xadmittances, theta, 0, 0, 0, 5.102E-2, 0.002));
		buses.add(new Bus(2, xadmittances, theta, 2, 0, 0, 1.502E-1, 0.033));
		buses.add(new Bus(1, xadmittances, theta, 4, -0.9, -0.5));
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = identifyNet(buses);

		double w0 = 1;
		double v0 = 1;
		double wi = 1;
		RealMatrix X0 = setUnknownValues(buses, deltaVoltageOrders, wi, buses.get(0).voltage.getValue());
		
		for (int i = 1; i < 4; i++) {
			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix(radmittances, xadmittances, wi);
			setActiveReactiveGen(buses, wi, w0, v0);
			
			X0 = setUnknownValues(buses, deltaVoltageOrders, wi,buses.get(0).voltage.getValue());
			
			ArrayList<ArrayList<DerivativeStructure>> pq = createEquations4(buses,
					Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances));
			
			ArrayList<Double> PQLossLoad = calculatePQLossLoad(cAdmittances, buses, wi);
			
			double[] mismatches = calculateMismatchMatrix(buses, wi, w0, v0, pq, PQLossLoad);

			RealMatrix fx0 = new Array2DRowRealMatrix(mismatches);
			System.out.println(fx0);
		
			double[][] Jacobian = constructJacabian2(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
					Admittance.createMadmittance(cAdmittances),Admittance.createTadmittance(cAdmittances),radmittances,xadmittances);
			
			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);
			
			for(int f=0;f<JJ.getRowDimension();f++) {
				for(int k=0;k<JJ.getColumnDimension();k++) {
					System.out.printf("%.2f\t",JJ.getEntry(f, k));
				}
				System.out.println();
			}
				
				
			RealMatrix X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
			
			updateUnknowns(X1, buses, deltaVoltageOrders, indexes, params, order);

			wi = X1.getEntry(X1.getRowDimension() - 2, 0);

			buses.get(0).voltage = new DerivativeStructure(params, order, 1, X1.getEntry(X1.getRowDimension() - 1, 0));
			
		
			
//			System.out.println("x0=\t"+X0);
//			System.out.println("fx0=\t"+fx0);
			
			for (int j = 0; j < X1.getRowDimension(); j++)
				System.out.printf("\t%s = %7.6f \t iteration %d %n", "Row".concat("" + j), X1.getEntry(j, 0), i);

			System.out.printf("%s%n", "--------------------------------------------");
			

		}

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
				b.p = 0.5*( ((1 / b.mp) * (w0 - wi)) + ((1 / b.nq) * (v0 - b.voltage.getValue())) );		
				b.q = 0.5*( ((1 / b.nq) * (v0 - b.voltage.getValue())) - ((1 / b.mp) * (w0 - wi)) ) ;
				System.out.printf("\tb.p=%f\t b.q=%f\t\n",b.p,b.q);
				System.out.printf("Voltage = %.2f\n",b.voltage.getValue());
				System.out.printf("Frequency = %.2f \n",wi);
			}
		}

	}

	public static ArrayList<Double> calculatePQLossLoad(Complex[][] cAdmittances, ArrayList<Bus> buses, double w) {
		ArrayList<Double> PQLossLoad = new ArrayList<>();

		
		double PLoad = 0.0;
		double QLoad = 0.0;
		double v0 = 1.01;
		double w0 = 1.0;
		double alpha = 1.0;
		double beta = 1.0;
		double kpf = 1.0;
		double kqf = -1.0;
		for (int k = 0; k < buses.size(); k++) {
			if (buses.get(k).type == 1) {
				PLoad += buses.get(k).p * Math.pow(buses.get(k).voltage.getValue() / v0, alpha) * (1 + kpf * (w - w0));
				QLoad += buses.get(k).q * Math.pow(buses.get(k).voltage.getValue() / v0, beta) * (1 + kqf * (w - w0));
			}
		}
		double PLoss = 0.0;
		double QLoss = 0.0;
		for (int k = 0; k < buses.size(); k++) {
				for (int n = 0; n < buses.size(); n++) {
					Complex temp0 = buses.get(k).cVolt.conjugate().multiply(buses.get(n).cVolt);
					Complex temp1 = temp0.add(buses.get(n).cVolt.conjugate().multiply(buses.get(k).cVolt));
					Complex temp2 = temp1.multiply(cAdmittances[k][n]);
//
//					 System.out.println("Volt K "+buses.get(k).cVolt);
//					 System.out.println("VoltageMultiplication :"+temp1);
//					 System.out.println("admittance :" +cAdmittances[k][n]);
//					 System.out.println("AdmittanceMultiplication :"+temp2);
//					 System.out.println();

					PLoss = PLoss + temp2.getReal();
					QLoss = QLoss + temp2.getImaginary();		
			}
		}
		
		PQLossLoad.add(0.5 * PLoss);
		PQLossLoad.add(-0.5 * QLoss);
		PQLossLoad.add(PLoad);
		
		PQLossLoad.add(QLoad);
		System.out.println("QLOSS ---->" + QLoss);
		System.out.println("QLoad ---->" + QLoad);
		System.out.println("PLOSS ---->" + PLoss);
		System.out.println("PLoad ---->" + PLoad);
		return PQLossLoad;
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

		System.out.printf("Ptotal:\t%f \nQtotal:\t %f\n", pTot, qTot);
		System.out.printf("Psys:\t%f \nQsys:\t %f\n", pSys, qSys);
		System.out.printf("Ptot-Psys:\t%f\n", pTot - pSys);
		System.out.printf("Qtot-Qsys:\t%f\n", qTot - qSys);
		System.out.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		// for(int i=0; i<mismatches.length;i++) {
		// System.out.print(mismatches[i]+"\t");
		// }
//		System.out.println();
		return mismatches;
	}
	
	public static double[][] constructJacabian2(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, ArrayList<Bus> buses, double w,
			Complex[][] cAdmittances, ArrayList<ArrayList<Integer>> indexes,RealMatrix admittance,RealMatrix theta,
			double [][]radmittances,double[][] xadmittances) {

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
				double admittanceR = -cAdmittances[k][n].pow(-1).getReal();;
				temp1 = -(FastMath.pow(admittanceX, 2) / w)
						/ (FastMath.pow((FastMath.pow(admittanceR, 2) + FastMath.pow(admittanceX, 2)), 3.0 / 2));
				double temp2 = 0;
				temp2 = -((admittanceX / (admittanceR * w)) / (1 + FastMath.pow(admittanceX / admittanceR, 2)));
				temp0 = temp0
						+ temp1 * buses.get(n).voltage.getValue()
								* FastMath.cos(buses.get(k).delta.getValue() - buses.get(n).delta.getValue()
										- theta.getEntry(k, n))
						+ temp2 * admittance.getEntry(k, n)
								* buses.get(n).voltage.getValue() * FastMath.sin(buses.get(k).delta.getValue()
										- buses.get(n).delta.getValue() - theta.getEntry(k, n));

			}

			j13[j] = temp0*buses.get(k).voltage.getValue();
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
						- temp2 * admittance.getEntry(k, n)
								* buses.get(n).voltage.getValue() * FastMath.cos(buses.get(k).delta.getValue()
										- buses.get(n).delta.getValue() - theta.getEntry(k, n));

			}

			j23[j] = temp0 * buses.get(k).voltage.getValue();
			j++;

		}

		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV] = j23[i - nofP];
		}
		double j33 = 0;
		for (Bus b : buses) {
			if (b.type == 2 && b.index!=0)
				j33 += -1 / (2*b.mp);
		}
		
		jacobian[nofP + nofQ][nofD + nofV] = j33;

		double j43 = 0;
		jacobian[nofP + nofQ + 1][nofP + nofQ] = j43;

		double j14[] = new double[nofP];
		j = 0;
		for (int k :pS) {
			j14[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k,0) * FastMath.cos(
					buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = 0; i < nofP; i++) {
			jacobian[i][nofD + nofV + 1] = j14[i];
		}

		double j24[] = new double[nofQ];
		j = 0;

		for (int k = 1; k < nofQ; k++) {

			j24[j] = buses.get(k).voltage.getValue() * admittance.getEntry(k, 0) * FastMath.sin(
					buses.get(k).delta.getValue() - buses.get(0).delta.getValue() - theta.getEntry(k, 0));
			j++;
		}
		for (int i = nofP; i < nofP + nofQ; i++) {
			jacobian[i][nofD + nofV + 1] = j24[i - nofP];
		}

		double j34 = 0;
		if(buses.get(0).type==2)
			jacobian[nofP + nofQ][nofD + nofV + 1] = -1 / (2*buses.get(0).nq);
		else
			jacobian[nofP + nofQ][nofD + nofV + 1] = j34;
		

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
				j42[j] = -1 / (2*buses.get(k).nq);
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
				if (buses.get(k).admittance.getEntry(k, i) != 0) {
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
			p.add(equationP.subtract(buses.get(k).p));
			if (buses.get(k).type != 0) {
				q.add(equationQ.subtract(buses.get(k).q));
			}
		}
		pq.add(p);
		pq.add(q);

		return pq;
	}

	
	public static RealMatrix setUnknownValues(ArrayList<Bus> buses,ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,double w,double v1){
		int nofD = deltaVoltageOrders.get(0).size();
		int nofV = deltaVoltageOrders.get(1).size();
		double[] deltas = new double[nofD];
		double[] voltages = new double[nofV];
		double[] unknowns = new double[nofD + nofV + 2];
		int i = 0;
		int j = 0;
		for (Bus b : buses) {
			if (b.index!=0) {
				if (b.type == 0) {
					deltas[i++] = b.delta.getValue();
				} else if (b.type == 1) {
					deltas[i++] = b.delta.getValue();
					voltages[j++] = b.voltage.getValue();
				} else if (b.type == 2) {
					deltas[i++] = b.delta.getValue();
					voltages[j++] = b.voltage.getValue();
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
			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders, ArrayList<ArrayList<Integer>> indexes, int params,
			int order) {
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
		
		for(Bus b: buses) {
			b.cVolt = Admittance.polarToComplex(b.voltage.getValue(), b.delta.getValue());
			
		}
		
	}

}
