package nRalgorithm;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;

public class Revision4 {

	public static void main(String[] args) {
		int nofB = 33;
		double[][] xadmittances = new double[nofB][nofB];
		double[][] radmittances = new double[nofB][nofB];
		File file = new File("case33.txt");
		try {
			Scanner scan = new Scanner(file);
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					xadmittances[i][j] = scan.nextDouble();

				}
			}
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					radmittances[i][j] = scan.nextDouble();
				}
			}
			scan.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;
		double w0 = 1.00;
		double v0 = 1.01;
		
		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;
		
		int params = 66;
		int order = 2;

		X1 = null;
		fx0 = null;
		wi = 1.0;

		prev_Mismatches = Double.MAX_VALUE;
		flag = true;
		RealMatrix prevX1 = null;
		RealMatrix prevfx0 = null;
		double dLoad = -0.03047;
		Jacobian = null;
		double[] position = {-0.764751,	-0.041600,	0.087224,	-0.044049,	
				-1.377096,	-3.257998,	-5.525477,	-1.790784};
		boolean flag3 = true;
		double vStart = 1;
		double dStart = -0.01;
		buses = new ArrayList<Bus>();

		buses.add(new Bus(1, 0, params,vStart,dStart,  dLoad));//1
		buses.add(new Bus(1, 2, params,vStart,dStart,  dLoad));//2
		buses.add(new Bus(1, 4, params,vStart,dStart,  dLoad));//3
		buses.add(new Bus(1, 6, params,vStart,dStart,  dLoad));//4
		buses.add(new Bus(1, 8, params,vStart,dStart,  dLoad));//5
		buses.add(new Bus(1, 10, params,vStart,dStart,  dLoad));//6
		buses.add(new Bus(1, 12, params, vStart,dStart, dLoad));//7
		buses.add(new Bus(1, 14, params,vStart,dStart,  dLoad));//8
		buses.add(new Bus(1, 16, params,vStart,dStart,  dLoad));//9
		buses.add(new Bus(1, 18, params,vStart,dStart,  dLoad));//10
		buses.add(new Bus(1, 20, params,vStart,dStart,  dLoad));//11
		buses.add(new Bus(1, 22, params,vStart,dStart,  dLoad));//12
		buses.add(new Bus(2, 24, params, 0 , 0, position[0],position[4], 1.2, 1.6,vStart,dStart)); //13
		buses.add(new Bus(1, 26, params, vStart,dStart, dLoad));//14
		buses.add(new Bus(2, 28, params, 0 , 0, position[1],position[5], 1.2, 1.6,vStart,dStart)); //15
		buses.add(new Bus(1, 30, params, vStart,dStart, dLoad));//16
		buses.add(new Bus(1, 32, params, vStart,dStart, dLoad));//17
		buses.add(new Bus(1, 34, params, vStart,dStart, dLoad));//18
		buses.add(new Bus(1, 36, params, vStart,dStart, dLoad));//19
		buses.add(new Bus(1, 38, params, vStart,dStart, dLoad));//19
		buses.add(new Bus(1, 40, params, vStart,dStart, dLoad));//21
		buses.add(new Bus(1, 42, params,vStart,dStart,  dLoad));//22
		buses.add(new Bus(1, 44, params, vStart,dStart, dLoad));//23
		buses.add(new Bus(1, 46, params,vStart,dStart,  dLoad));//24
		buses.add(new Bus(2, 48, params, 0 , 0, position[2],position[6], 1.2, 1.6,vStart,dStart)); //25
		buses.add(new Bus(1, 50, params, vStart,dStart, dLoad));//26
		buses.add(new Bus(1, 52, params,vStart,dStart,  dLoad));//27
		buses.add(new Bus(1, 54, params, vStart,dStart, dLoad));//28
		buses.add(new Bus(1, 56, params, vStart,dStart, dLoad)); //29
		buses.add(new Bus(1, 58, params, vStart,dStart, dLoad));//30
		buses.add(new Bus(1, 60, params, vStart,dStart, dLoad));//31
		buses.add(new Bus(1, 62, params, vStart,dStart, dLoad));//32
		buses.add(new Bus(2, 64, params, 0 , 0, position[3],position[7], 1.2, 1.6,vStart,dStart)); //33
		
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);
		
		innerloop: for (int i = 1; i <= 150 && flag; i++) {  


			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

			ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

			PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

			RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
					buses.get(0).voltage.getValue());

			RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

			RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
					tAdmittances);

			double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

			fx0 = new Array2DRowRealMatrix(mismatches);

			Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
					Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
					radmittances, xadmittances);

			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

			try {
				try {

					X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
					if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
						System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
						break innerloop;
					}
				} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

					e.printStackTrace();
				}
			} catch (SingularMatrixException e) {

				for (int f = 0; f < Jacobian.length; f++) {
					for (int j = 0; j < Jacobian[1].length; j++) {
						System.out.printf("%.2f \t", Jacobian[f][j]);
					}
					System.out.println();
				}
				for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
					for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
						System.out.printf("%.2f \t", Admittance.createMadmittance(cAdmittances).getEntry(f, j));
					}
					System.out.println();
				}

				for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
					for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
						System.out.printf("%.2f,%.2f \t", cAdmittances[f][j].getReal(),
								cAdmittances[f][j].getImaginary());
					}
					System.out.println();
				}

				break innerloop;

			}
			if (ModifiedNR.sumMatrix(fx0) <= 1E-4) {

				System.out.println("\nFull convergence achieved. Exiting..." + "*".repeat(10));
				flag = false;
				System.out.println(X1);

				
			} else if (prev_Mismatches + 0.5e10 <= ModifiedNR.sumMatrix(fx0)) {

				System.out.printf("\nMissed the local optimum. Exiting... \t At iteration : %d", i);
				System.out.printf("\nRejected sum of mismatches: %.5f", ModifiedNR.sumMatrix(fx0));
				X1 = prevX1;
				flag = false;
			} else {

				prevX1 = X1;
				prevfx0 = fx0;
				prev_Mismatches = ModifiedNR.sumMatrix(fx0);

				wi = X1.getEntry(X1.getRowDimension() - 2, 0);

				buses.get(0).voltage = new DerivativeStructure(params, order, 1,
						X1.getEntry(X1.getRowDimension() - 1, 0));

				ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
			}

			System.out.printf("%.5f,",ModifiedNR.sumMatrix(fx0));
			
			
		}	
	}

}


