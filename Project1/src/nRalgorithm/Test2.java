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

public class Test2 {

	public static void main(String[] args) {

		int nofB = 6;
		double[][] xadmittances = new double[nofB][nofB];
		double[][] theta = new double[nofB][nofB];
		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		ArrayList<Double> PQLossLoad = null;
		double treat = 1.0;

		File file = new File("/Users/my_mac/git/com.yavuz.git.first/Project1/test.txt");
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
					theta[i][j] = Admittance.getAng(new Complex(radmittances[i][j], xadmittances[i][j]));
				}
			}
			scan.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		RealMatrix X1 = null;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		// PV BUS
		buses.add(new Bus(1, 0, theta.length * 2, -0.03525369 * treat, -0.026188417 * treat));
//		buses.add(new Bus(1, 0, theta.length * 2, 0, 0, 1.0));
		// Droop Bus 9.4E-4 0.0013
		buses.add(new Bus(2, 2, theta.length * 2, 0, 0, 0.48794, 0.98998, 0.9));
		buses.add(new Bus(2, 4, theta.length * 2, 0, 0, 0.41974, 0.69433, 0.9));
		buses.add(new Bus(2, 6, theta.length * 2, 0, 0, 0.82840, 0.99666, 0.9));
		
//		buses.add(new Bus(2, 2, theta.length * 2, 0, 0, 9.4E-3, 0.0013));
//		buses.add(new Bus(2, 4, theta.length * 2, 0, 0, 9.4E-3, 0.0013));
//		buses.add(new Bus(2, 6, theta.length * 2, 0, 0, 9.4E-3, 0.0013));

		// PQ Bus
//		buses.add(new Bus(1, 8, theta.length * 2, -0.03525369 * treat, -0.06188417 * treat));
		buses.add(new Bus(1, 8, theta.length * 2, -0.0 * treat, -0.0 * treat));
		buses.add(new Bus(1, 10, theta.length * 2, -0.04417614 * treat, -0.008281924 * treat));
//		buses.add(new Bus(1, 12, theta.length * 2, -0.4417614 * treat, -0.08281924 * treat));
		//PV Bus
//		buses.add(new Bus(0, 12, theta.length * 2, 0.02, 0.000, 1.01));
		
//		buses.add(new Bus(2, 14, theta.length * 2,0, 0,.09492, .06811, 0.9));
		
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

		int params = buses.size() * 2;
		int order = 2;
		RealMatrix fx0 = null;
		double[][] Jacobian = null;
		double wi = 1.0;

		double prev_Mismatches = Double.MAX_VALUE;
		RealMatrix prevX1 = new Array2DRowRealMatrix(12, 1);
		boolean flag = true;
		double tolerance = 0.02;
		RealMatrix prevfx0 = null;

		long cur = System.currentTimeMillis();

		innerloop: for (int i = 1; i <= 400 && flag; i++) {
			double w0 = 1.0;
			double v0 = 1.01;

			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

			ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

			RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi, buses.get(0).voltage.getValue());

			PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses,
					Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances));

			double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

			fx0 = new Array2DRowRealMatrix(mismatches);

			Jacobian = ModifiedNR.constructJacabian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
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

			if (prev_Mismatches + tolerance < ModifiedNR.sumMatrix(fx0)) {
				System.out.println("Missed the local optimum. Exiting...");
				System.out.printf("Rejected sum of mismatches: %.5f", ModifiedNR.sumMatrix(fx0));
				X1 = prevX1;
				flag = false;

			} else if (ModifiedNR.sumMatrix(fx0) < 1E-5) {
				System.out.println("Full convergence achieved. Exiting...");
				flag = false;
			} else {
				System.out.println("Difference: " + X1.subtract(prevX1));
				prevX1 = X1;
				prevfx0 = fx0;
				prev_Mismatches = ModifiedNR.sumMatrix(fx0);

				wi = X1.getEntry(X1.getRowDimension() - 2, 0);

				buses.get(0).voltage = new DerivativeStructure(params, order, 1,
						X1.getEntry(X1.getRowDimension() - 1, 0));

				ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, indexes, params, order);
			}

			System.out.printf("\nLatest sum of mismatches = %.5f \n\n", prev_Mismatches);
			for (int j = 0; j < X1.getRowDimension(); j++) {
				if (j < deltaVoltageOrders.get(0).size())
					System.out.printf("\t%s = %7.6f \t iteration %d %n", "Row".concat("" + j),
							X1.getEntry(j, 0) * 180 / Math.PI, i);
				else
					System.out.printf("\t%s = %7.6f \t iteration %d %n", "Row".concat("" + j), X1.getEntry(j, 0), i);
			}
			System.out.printf("%s%n", "------------------------------------------------");

		}
		if (!Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
			System.out.println("Total Time Elapsed (in msec) : " + (System.currentTimeMillis() - cur));
			for (Bus b : buses) {
				System.out.printf(
						"\nBus index: %d \t Bus type: %s\n" + "Bus Voltage: %.4f\n" + "Bus Angle: %.4f\n"
								+ "Bus Active Power: %.4f\n" + "Bus Reactive Power: %.4f\n"
								+ "-------------------------------------------------",
						buses.indexOf(b), b.type == 0 ? "PV" : b.type == 1 ? "PQ" : "DROOP", b.voltage.getValue(),
						b.delta.getValue()* 180 / Math.PI , b.p, b.q);

			}

			for (Bus b : buses)
				if (b.type == 2)
					System.out.printf("\nDroop coefficients of bus %d  \t mp: %.5f nq: %.5f\n", buses.indexOf(b), b.mp,
							b.nq);

			System.out.println("\nFinal Mismatches = " + prevfx0);
			System.out.println("\nMismatch Summation: " + prev_Mismatches);
			System.out.printf("\nPLoss: %.5f\tPLoad: %.5f\n" + "QLoss: %.5f\tQLoad: %.5f \n", -PQLossLoad.get(0),
					PQLossLoad.get(2), -PQLossLoad.get(1), PQLossLoad.get(3));
			System.out.printf("\nTotal Active Power Demand: %.5f\n" + "Total Reactive Power Demand: %.5f \n",
					-PQLossLoad.get(0) + PQLossLoad.get(2), -PQLossLoad.get(1) + PQLossLoad.get(3));
			System.out.printf("\nSteady state system frequency: %.5f\n", wi);
			System.out.printf("Steady state reference voltage value: %.5f\n", buses.get(0).voltage.getValue());

		}

	}

}
