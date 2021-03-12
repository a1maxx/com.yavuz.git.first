package nRalgorithm;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
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

public class ModifiedNR6 {

	public static void main(String[] args) {
		double[][] xadmittances = new double[3][3];
		double[][] theta = new double[3][3];
		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		ArrayList<Double> PQLossLoad = null;

		File file = new File("test4.txt");
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
		RealMatrix fx0 = null;

		double wi = 1.0;
		RealMatrix prevX1 = new Array2DRowRealMatrix(6, 1);
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		
		double[][] Jacobian = null;
		boolean flag2 = true;
		int m = 0;
		int u = 0;
		ArrayList<Double> initialD = new ArrayList<Double>();
		ArrayList<Double> initialV = new ArrayList<Double>();
		long cur = System.currentTimeMillis();
		while (u < 30 && !((System.currentTimeMillis()-cur)>1.5e4 && u<3) )  {

			int params = 6;
			int order = 2;
			X1 = null;
			fx0 = null;
			wi = 1.0;
			prevfx0 = null;
			prevX1 = new Array2DRowRealMatrix(6, 1);
			prev_Mismatches = Double.MAX_VALUE;
			flag = true;
			Jacobian = null;

			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, 6, -0.5, -0.2));
			buses.get(0).delta = new DerivativeStructure(6, 2, 0, 0.0);

			buses.add(new Bus(2, 2, 6, 0, 0, 0.05, 0.85, 3.0));
			buses.add(new Bus(2, 4, 6, 0, 0, 0.09, 0.28, 1.5));

			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			int ITERATION = 300;

			innerloop: for (int i = 1; i <= ITERATION && flag; i++) {
				double w0 = 1.00;
				double v0 = 1.01;

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
				if (ModifiedNR.sumMatrix(fx0) < 1E-6) {
					boolean flag3 = true;

					for (Bus b : buses) {
						if (b.type == 2 && (b.p < 0 || b.q < 0)) {
							flag3 = false;
						}
					}
					if (flag3) {
						System.out.println("\nFull convergence achieved. Exiting...");

						for (Bus b : buses) {
							initialD.add(b.initialValue_d);
							initialV.add(b.initialValue_v);
						}

						flag = false;
						flag2 = false;
						u++;

					}

				} else if (prev_Mismatches + 10.00984  <= ModifiedNR.sumMatrix(fx0)) {
//
//					System.out.printf("\nMissed the local optimum. Exiting... \t At iteration : %d", i);
//					System.out.printf("\nRejected sum of mismatches: %.5f", ModifiedNR.sumMatrix(fx0));

					m++;

					if (m % 100 == 0) {
						System.out.println("");
					} else {
						System.out.print(u);
					}

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

				if (i == ITERATION) {
					System.out.print("<I>");
				}

			}
		}

		int f = 0;
		outerloop: for (int i = 0; i < initialD.size(); i++) {
			for (int j = 0; j < 3; j++) {
				System.out.printf("%.5f\t", initialD.get(f));
				if (f + 1 == initialD.size())
					break outerloop;
				else
					f++;
			}
			System.out.println();
		}
		
		System.out.println("\n"+"-".repeat(100));
		f = 0;
		
		outerloop: for (int i = 0; i < initialV.size(); i++) {
			for (int j = 0; j < 3; j++) {
				System.out.printf("%.5f\t", initialV.get(f));
				if (f + 1 == initialV.size())
					break outerloop;
				else
					f++;

			}
			System.out.println();
		}

//		if (!Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
//			System.out.println("\nTotal Time Elapsed (in msec) : " + (System.currentTimeMillis() - cur));
//			for (Bus b : buses) {
//				System.out.printf(
//						"\nBus index: %d \t Bus type: %s\n" + "Bus Voltage: %.4f\n" + "Bus Angle: %.4f\n"
//								+ "Bus Active Power: %.4f\n" + "Bus Reactive Power: %.4f\n"
//								+ "-------------------------------------------------",
//						buses.indexOf(b), b.type == 0 ? "PV" : b.type == 1 ? "PQ" : "DROOP", b.voltage.getValue(),
//						b.delta.getValue(), b.p, b.q);
//
//			}
//
//			for (Bus b : buses)
//				if (b.type == 2)
//					System.out.printf("\nDroop coefficients of bus %d  \t mp: %.5f nq: %.5f\n", buses.indexOf(b), b.mp,
//							b.nq);
//
//			System.out.println("\nFinal Mismatches = " + prevfx0);
//			System.out.println("\nMismatch Summation: " + prev_Mismatches);
//			System.out.printf("\nPLoss: %.5f\tPLoad: %.5f\n" + "QLoss: %.5f\tQLoad: %.5f \n", -PQLossLoad.get(0),
//					PQLossLoad.get(2), -PQLossLoad.get(1), PQLossLoad.get(3));
//			System.out.printf("\nTotal Active Power Demand: %.5f\n" + "Total Reactive Power Demand: %.5f \n",
//					-PQLossLoad.get(0) + PQLossLoad.get(2), -PQLossLoad.get(1) + PQLossLoad.get(3));
//			System.out.printf("\nSteady state system frequency: %.5f\n", wi);
//			System.out.printf("Steady state reference voltage value: %.5f\n", buses.get(0).voltage.getValue());
//		}

	}

}
