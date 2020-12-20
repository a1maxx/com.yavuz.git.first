package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class Test2 {

	public static void main(String[] args) {

		// Ybus elements
		double Y11, Y12, Y13, Y21, Y22, Y23, Y31, Y32, Y33;
		Y11 = 14;
		Y12 = 10;
		Y13 = 4;
		Y21 = 10;
		Y22 = 15;
		Y23 = 5;
		Y31 = 4;
		Y32 = 5;
		Y33 = 9;

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
		double[][] admittances = { { Y11, Y12, Y13 }, { Y21, Y22, Y23 }, { Y31, Y32, Y33 } };
		double[][] theta = { { t11, t12, t13 }, { t21, t22, t23 }, { t31, t32, t33 } };

		ArrayList<Bus> buses = new ArrayList<Bus>();
		buses.add(new Bus(3, admittances, theta, 0, 0, 0));
		buses.add(new Bus(0, admittances, theta, 2, -0.9, -0.5));
		buses.add(new Bus(0, admittances, theta, 4, 0.6, 0));
		buses.get(0).voltage = new DerivativeStructure(params, order, 1, 1.0);
		buses.get(0).delta = new DerivativeStructure(params, order, 0, 0);
		buses.get(2).voltage = new DerivativeStructure(params, order, 5, 1.01);
		
		
		// Creating equations checked correct...
//			equations = Bus.createEquations2(buses);
//			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations3(buses);
//			for(int i =0;i<equations.length;i++)
//				System.out.println(equations[i].getValue());
//			for(int j=0;j<2;j++)
//				for(int k=0;k<pq.get(j).size();k++)
//			System.out.println(pq.get(j).get(k).getValue());
		
		// Creating orders checked correct...
//			ArrayList<Integer[]> orders = GenericNR.createOrders(buses);
//			for(Integer [] a:orders) {
//				for(int i =0;i<a.length;i++) {
//					System.out.print(a[i]);
//				}
//				System.out.println("-------");
//			}
//			buses.get(1).type=1;
//			buses.get(2).type=0;
		
//			
//			for(int k=0 ;k<2;k++) {
//				ArrayList<Integer[]> orders3 = orders2.get(k);
//			for(Integer [] a:orders3) {
//				for(int i =0;i<a.length;i++) {
//					System.out.print(a[i]);
//				}
//				System.out.println();
//			}
//			}
		
//		 //Setting Unknowns checked Correct...
//			ArrayList<Double> unknowns0 = GenericNR.calculateUnknowns(buses);
//			for(Double d:unknowns0) {
//				System.out.println(d);
//			}
//			buses.get(1).type=1;
//			buses.get(2).type=0;
//			ArrayList<ArrayList<Integer[]>> orders2 = ModifiedNR.createOrders2(buses);
//			RealMatrix x0 = ModifiedNR.setUnknownValues(buses, orders2, 1.0, buses.get(0).voltage.getValue());
//			System.out.println(x0);
		
		
			
//			X0 = new Array2DRowRealMatrix(GenericNR.convertArray(unknowns0));
			
//			double[][] jacobs = GenericNR.calculateJacobian(orders, equations);
//			RealMatrix J = new Array2DRowRealMatrix(jacobs);
//			double[] functions0 = GenericNR.calculateFunctions(equations);
//			RealMatrix fx0 = new Array2DRowRealMatrix(functions0);
//			RealMatrix JI = MatrixUtils.inverse(J);
//			X0 = X0.subtract(JI.multiply(fx0));
//			GenericNR.updateUnknowns(X0, params, order, orders, buses);
			
			DerivativeStructure[] equations = Bus.createEquations2(buses);
			ArrayList<Integer[]> orders = GenericNR.createOrders(buses);
			double[][] jacobs = GenericNR.calculateJacobian(orders, equations);
			RealMatrix J = new Array2DRowRealMatrix(jacobs);
			System.out.println(J);
			
			//System.out.println(buses.get(0).admittance);
			
			buses.get(1).type=1;
			buses.get(2).type=1;
			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);
			
			
			
	double X11, X12, X13, X21, X22, X23, X31, X32, X33;
	X11 = 0;
	X12 = 0.1;
	X13 = 0.25;
	X21 = 0.1;
	X22 = 0;
	X23 = 0.2;
	X31 = 0.25;
	X32 = 0.2;
	X33 = 0;
	double wi=1.0;
	int N=3;
			double[][] xadmittances = { { X11, X12, X13}, { X21, X22, X23 }, { X31, X32, X33 } };
			double [][] radmittances = new double[xadmittances.length][xadmittances[1].length];
			for (int i = 0, len = radmittances.length; i < len; i++)
			    Arrays.fill(radmittances[i], 0);
			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix(radmittances, xadmittances,wi);
			ArrayList<ArrayList<DerivativeStructure>> pq2 = ModifiedNR.createEquations4(buses,Admittance.createMadmittance(cAdmittances),
					Admittance.createTadmittance(cAdmittances));
			
			double [][] mAdmittance = new double[buses.size()][buses.size()];
			double [][] tAdmittance = new double[buses.size()][buses.size()];
			
			for(int i=0;i<cAdmittances.length;i++) {
				for(int j=0;j<cAdmittances[0].length;j++) {
					mAdmittance[i][j]=Admittance.getMag(cAdmittances[i][j]);
					tAdmittance[i][j] =Admittance.getAng(cAdmittances[i][j]);
				}
			}
			
			
			// Jacobian creation checked correct...
			double [][] Jacobian = ModifiedNR.constructJacabian(deltaVoltageOrders, pq2, buses, wi, cAdmittances, indexes);
			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);
			int [] r= {0, 1, 2,3};
			int [] c = {0, 1, 2,3};
			
			System.out.println(pq2.get(1).get(0).getPartialDerivative(0,0,0,0,0,1));
			System.out.println(equations[1].getPartialDerivative(0,0,0,0,0,1));
			
			System.out.println(pq2.get(1).get(0).getPartialDerivative(0,0,0,1,0,0));
			System.out.println(pq2.get(1).get(0).getPartialDerivative(ArrayUtils.toPrimitive(deltaVoltageOrders.get(1).get(0))));
			for(int i:deltaVoltageOrders.get(0).get(1)) {
				System.out.printf("%d\t",i);
			}
			System.out.println();
			System.out.println(equations[2].getPartialDerivative(0,0,0,1,0,0));
			
			System.out.println(JJ.getSubMatrix(r,c));
			//System.out.println(MatrixUtils.inverse(JJ));
			System.out.println(MatrixUtils.inverse(J));
			
			
			
			

	}

		// TODO Auto-generated method stub

	}


