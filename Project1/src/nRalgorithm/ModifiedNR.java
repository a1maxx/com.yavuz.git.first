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

		
		double[][] admittances = { { Y11, Y12, Y13}, { Y21, Y22, Y23 }, { Y31, Y32, Y33 } };
		double[][] theta = { { t11, t12, t13 }, { t21, t22, t23 }, { t31, t32, t33 } };
		
		
		ArrayList<ArrayList<Admittance>> admittanceS = new ArrayList<>();
		double [][] xadmittances = new double[admittances.length][admittances[1].length];
		for (int i = 0, len = xadmittances.length; i < len; i++)
		    Arrays.fill(xadmittances[i], 0.002);
		double w = 1.0;
		Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix(admittances, xadmittances,w);
		
		for (int i = 0; i < 3; i++) {
			ArrayList<Admittance> temp0 = new ArrayList<>();
			for (int j = 0; j < 3; j++) {
				Admittance temp1 = new Admittance(0.002, 0.0001);
				temp0.add(temp1);
			}
			admittanceS.add(temp0);
		}
		
		ArrayList<Bus> buses = new ArrayList<Bus>();
		buses.add(new Bus(2, admittances, theta, 0, 0, 0,1.02,0));
		buses.add(new Bus(2, admittances, theta, 2, 0, 0,1.01,0));
		buses.add(new Bus(1, admittances, theta, 4, -0.6, -0.372,0.97,0.0014));
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes =identifyNet(buses);
		double w0 = 1;
		double v0 = 1;
		double wi = 1;
		double v1 = 1;
		RealMatrix X0 = setUnknownValues(buses, deltaVoltageOrders, wi, v1);
		int N =3;
		for(int i = 1;i<10;i++) {
			setActiveReactiveGen(buses,wi,w0,v0);
			cAdmittances = Admittance.constructComplexAdmittanceMatrix(admittances, xadmittances,wi);
			ArrayList<ArrayList<DerivativeStructure>> pq2 = createEquations3(buses);
			ArrayList<Double> PQLossLoad2 = calculatePQLossLoad(cAdmittances,buses);
			double [] mismatches2 = calculateMismatchMatrix(buses, wi, w0, v0, pq2, PQLossLoad2);
			RealMatrix fx0 = new Array2DRowRealMatrix(mismatches2);
		
			double [][] Jacobian = constructJacabian(deltaVoltageOrders, pq2, N, buses, wi, admittanceS, indexes);
			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);
			RealMatrix X1= X0.add(MatrixUtils.inverse(JJ).multiply(fx0));
			updateUnknowns(X1,buses,deltaVoltageOrders,indexes,params,order);
			//Update of w and V
			wi = X1.getEntry(X1.getRowDimension()-2, 0);
			buses.get(0).voltage = new DerivativeStructure(params, order, 0, X1.getEntry(X1.getRowDimension()-1, 0));
			
			for(int j=0;j<X1.getRowDimension();j++) 
				System.out.printf("\t%s = %7.6f \t iteration %d %n", "Column".concat(""+j) ,X1.getEntry(j, 0), i);
			
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

	public static void setActiveReactiveGen(ArrayList<Bus> buses, double wi, double w0, double v0) {
		for (Bus b : buses) {
			if (b.type == 2) {
				b.p = 1 / b.mp * (w0 - wi);
				b.q = 1 / b.nq * (v0 - b.voltage.getValue());
			}
		}

	}

	public static double[] calculateMismatchMatrix(ArrayList<Bus> Buses, double w, double w0, double v0,
			ArrayList<ArrayList<DerivativeStructure>> pq,ArrayList<Double> PQLossLoad) {

		double pSys = 0;
		double qSys = 0;
		double pTot = PQLossLoad.get(0)+PQLossLoad.get(2);
		double qTot = PQLossLoad.get(1)+PQLossLoad.get(3);
		for (Bus b : Buses) {
			if (b.type == 2) {
				pSys += (1 / b.mp) * (w0 - w);
				qSys += (1 / b.nq) * (v0 - b.voltage.getValue());
			}
		}
		
		int nofP =pq.get(0).size();
		int nofQ =pq.get(1).size();
		
		
		double[] pMismatch = new double[nofP];
		double[] qMismatch = new double[nofQ];
		for (int i = 0; i < pq.get(0).size(); i++) 
				pMismatch[i] = pq.get(0).get(i).getValue();
		
		for (int i = 0; i < pq.get(1).size(); i++) 
			qMismatch[i] = pq.get(1).get(i).getValue();		
		
		
		double[] mismatches = new double[nofP+nofQ+2];
		System.arraycopy(pMismatch, 0, mismatches, 0, pMismatch.length);
		System.arraycopy(qMismatch, 0, mismatches, pMismatch.length, qMismatch.length);
		mismatches[(pMismatch.length + qMismatch.length)] = pTot - pSys;
		mismatches[(pMismatch.length + qMismatch.length) + 1] = qTot - qSys;
		
		System.out.printf("Ptotal:\t%f \nQtotal:\t %f\n",pTot,qTot);
		System.out.printf("Psys:\t%f \nQsys:\t %f\n",pSys,qSys);
		System.out.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		return mismatches;
	}

	public static double[][] constructJacabian(ArrayList<ArrayList<Integer[]>> deltaVoltageOrders,
			ArrayList<ArrayList<DerivativeStructure>> pq, int N, ArrayList<Bus> buses, double w,
			ArrayList<ArrayList<Admittance>> admittance,ArrayList<ArrayList<Integer>> indexes) {
		
		int nofP =pq.get(0).size();
		int nofQ =pq.get(1).size();
		int nofD =deltaVoltageOrders.get(0).size();
		int nofV = deltaVoltageOrders.get(1).size();
//		N = nofP+ nofQ;
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
		
		double j23[] = new double[nofQ];
		ArrayList<Integer> qS = new ArrayList<Integer>();
		qS.addAll(indexes.get(1));
		qS.addAll(indexes.get(2));
		j = 0;
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
		if (buses.get(0).type == 2)
			j44 = -1 / buses.get(0).nq;
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
		
	public static ArrayList<Double> calculatePQLossLoad(Complex[][] cAdmittances,ArrayList<Bus> buses) {
		ArrayList<Double> PQLossLoad = new ArrayList<>();
		
		double PLoss=0;
		double QLoss=0;
		double PLoad=0;
		double QLoad=0;
			for(int k=0;k<buses.size();k++) {
				if(buses.get(k).type==1){
					PLoad+=buses.get(k).p;
					QLoad+=buses.get(k).q;
				}
				for(int n=0;n<buses.size();n++) {
					Complex temp0= buses.get(k).cVolt.conjugate().multiply(buses.get(n).cVolt);
					Complex temp1= temp0.add(buses.get(n).cVolt.conjugate().multiply(buses.get(k).cVolt));
					Complex temp2 =temp1.multiply(cAdmittances[k][n]);
					PLoss = PLoss + temp2.getReal();
					QLoss = QLoss + temp2.getImaginary();
				}
			}
		
		PQLossLoad.add(0.5*PLoss);
		PQLossLoad.add(-0.5*QLoss);
		PQLossLoad.add(PLoad);
		PQLossLoad.add(QLoad);
		return PQLossLoad;
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
		
		//for(Bus b: buses) {
		//	b.cVolt = new Complex(b.voltage.getValue(),b.delta.getValue());
		//}
		
	}

}
