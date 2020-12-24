package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class ModifiedNR2 {

	public static void main(String[] args) {

		double X11, X12, X13, X21, X22, X23, X31, X32, X33;
		X11 = 0;
		X12 = 0.1;
		X13 = 1.25;
		X21 = 0.1;
		X22 = 0;
		X23 = 0.2;
		X31 = 1.25;
		X32 = 0.2;
		X33 = 0;
		
		double X14,X24,X34,X44,X41,X42,X43;
		X14=0.2;
		X24 =0.1;
		X34 =0.7;
		X44= 0;
		X43 =0.7;
		X42=0.1;
		X41=0.1;
		

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

		double t14, t24, t34, t44, t41, t42, t43;
		t14 = pi / 2;
		t24 = pi / 2;
		t34 = pi / 2;
		t44 = -pi / 2;
		t43 = pi / 2;
		t42 = pi / 2;
		t41 = pi / 2;


		// Known variables

	
		double[][] xadmittances = { { X11, X12, X13 ,X14}, { X21, X22, X23,X24 }, { X31, X32, X33 ,X34 },{X41,X42,X43,X44} };
		double[][] theta = { { t11, t12, t13 ,t14}, { t21, t22, t23,t24 }, { t31, t32, t33 ,t34 },{t41,t42,t43,t44} };

		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		for (int i = 0, len = radmittances.length; i < len; i++)
			Arrays.fill(radmittances[i], 0);

		ArrayList<Bus> buses = new ArrayList<Bus>();
		//PV BUS
		buses.add(new Bus(0,0,theta.length*2,0.2,0.1,1.0));
		//Droop Bus
		buses.add(new Bus(2,2,theta.length*2,0, 0, 0.0692791530265256, 0.40898603960680735));
		//PQ Bus
		buses.add(new Bus(1,4,theta.length*2, -0.8, -0.1));	
		buses.add(new Bus(1, 6,theta.length*2, -0.1, -0.1));
		
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);
		
		int params = buses.size()*2;
		int order = 2;
	
		double wi = 1;
		long cur =System.currentTimeMillis();
		for (int i = 1; i < 10; i++) {
			double w0 = 1.0;
			double v0 = 1.1;
			
			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix(radmittances, xadmittances, wi);
			
			ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);
			
			RealMatrix X0  = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,buses.get(0).voltage.getValue());
			
			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations4(buses,
					Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances));
			
			ArrayList<Double> PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);
			
			double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, wi, w0, v0, pq, PQLossLoad,indexes);

			RealMatrix fx0 = new Array2DRowRealMatrix(mismatches);
			
			double[][] Jacobian = ModifiedNR.constructJacabian2(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
					Admittance.createMadmittance(cAdmittances),Admittance.createTadmittance(cAdmittances),radmittances,xadmittances);
			
			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);
			
//			for(int f=0;f<JJ.getRowDimension();f++) {
//				for(int k=0;k<JJ.getColumnDimension();k++) {
//					System.out.printf("%.2f\t",JJ.getEntry(f, k));
//				}
//				System.out.println();
//			}	
			
			RealMatrix X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
//			System.out.println(X1);
			

			wi = X1.getEntry(X1.getRowDimension() - 2, 0);

			buses.get(0).voltage = new DerivativeStructure(params, order, 1, X1.getEntry(X1.getRowDimension()-1, 0));
			
			ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, indexes, params, order);
		
			System.out.println("X0=\t"+X0);
			System.out.println("X1=\t"+X1);
			System.out.println("fx0=\t"+fx0);
			
			for (int j = 0; j < X1.getRowDimension(); j++)
				System.out.printf("\t%s = %7.6f \t iteration %d %n", "Row".concat("" + j), X1.getEntry(j, 0), i);

			System.out.printf("%s%n", "--------------------------------------------");
			

		}
		System.out.println("Total Time Elapsed : " + (System.currentTimeMillis()-cur));
	}

}
