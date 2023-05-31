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

public class Revision1 {

	public static void main(String[] args) {
		int nofB = 6;
		double[][] xadmittances = new double[nofB][nofB];
		double[][] radmittances = new double[nofB][nofB];
		File file = new File("case6ww.txt");
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
		
		int params = 12;
		int order = 2;

		X1 = null;
		fx0 = null;
		wi = 1.0;

		prev_Mismatches = Double.MAX_VALUE;
		flag = true;
		RealMatrix prevX1 = null;
		RealMatrix prevfx0 = null;

		Jacobian = null;

		
		boolean flag3 = true;
		double vStart = 0.5;
		double dStart = 0.1;
		buses = new ArrayList<Bus>();
		buses.add(new Bus(1, 0, params,vStart,dStart, -0.01 ));//1
		buses.add(new Bus(1, 2, params,vStart,dStart,  -0.01));//2
		buses.add(new Bus(1, 4, params,vStart,dStart, -0.01));//3
		buses.add(new Bus(2,6,params,0 ,0 , 0.5,0.5,1.2,1.6,vStart,dStart)); //4
		buses.add(new Bus(2, 8, params,0 , 0, 0.5, 0.5, 1.2, 1.6,vStart,dStart)); //5
		buses.add(new Bus(2, 10, params, 0 ,0, 0.5, 0.5, 1.2, 1.6,vStart,dStart)); //6

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);
		
		innerloop: for (int i = 1; i <= 100 && flag; i++) {  


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

					RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
					double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
							radmittances, xadmittances, DELTA, X0);

					RealMatrix aj = new Array2DRowRealMatrix(arJacob);
					DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
					X1 = X0.add(DELTA);

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

				
			} else if (prev_Mismatches + 0.5e1 <= ModifiedNR.sumMatrix(fx0)) {

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

		}	
	}

}
