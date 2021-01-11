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

public class ModifiedNR2 {

	public static void main(String[] args) {

		Random random = new Random();

		double[][] xadmittances = new double[6][6];
		double[][] theta = new double[6][6];
		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		ArrayList<Double> PQLossLoad = null;

		File file = new File("/Users/my_mac/git/com.yavuz.git.first/Project1/test2.txt");
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
		double [] best_mp = new double[3];
		double [] best_nq = new double[3];
 		double min_mismatch = Double.MAX_VALUE;
 
 		int n_root = 3000;
		double[] mp = new double[n_root];
		double[] nq = new double[n_root];
		for (int i = 0; i < mp.length; i++)
			mp[i] = random.nextDouble();
		for (int i = 0; i < nq.length; i++)
			nq[i] = random.nextDouble();
//		outerloop:
		for (int m = 0; m < mp.length; m++) {
			
			RealMatrix X1 = null;
			ArrayList<Bus> buses = new ArrayList<Bus>();
			// PV Buses
			buses.add(new Bus(0, 0, theta.length * 2, 0, 0, 1.0));
			// Droop Buses
			buses.add(new Bus(2, 2, theta.length * 2, 0, 0, mp[random.nextInt(mp.length)], nq[random.nextInt(nq.length)]));
			buses.add(new Bus(2, 4, theta.length * 2, 0, 0, mp[random.nextInt(mp.length)], nq[random.nextInt(nq.length)]));
			buses.add(new Bus(2, 6, theta.length * 2, 0, 0, mp[random.nextInt(mp.length)], nq[random.nextInt(nq.length)]));
			// PQ Buses
			buses.add(new Bus(1, 8, theta.length * 2, -0.0825, -0.0618));
			buses.add(new Bus(1, 10, theta.length * 2, -0.0244, -0.0828));


			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			int params = buses.size() * 2;
			int order = 2;
			RealMatrix fx0 = null;
			double[][] Jacobian = null;
			double wi = 1;
			long cur = System.currentTimeMillis();
			innerloop:
			for (int i = 1; i <= 75; i++) {
				double w0 = 1.0;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

				RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

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
						if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))
								|| Math.abs(X1.getEntry(X1.getRowDimension() - 2, 0))
										+ Math.abs(X1.getEntry(X1.getRowDimension() - 1, 0)) > 3
								|| ModifiedNR.sumMatrix(fx0) >  1) {
							System.out.printf("\nDIVERGED at iteration %d", m);
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

				wi = X1.getEntry(X1.getRowDimension() - 2, 0);

				buses.get(0).voltage = new DerivativeStructure(params, order, 1,
						X1.getEntry(X1.getRowDimension() - 1, 0));

				ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, indexes, params, order);

			}
		
			if (!(Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))
					|| Math.abs(X1.getEntry(X1.getRowDimension() - 2, 0))
					+ Math.abs(X1.getEntry(X1.getRowDimension() - 1, 0)) > 3
					|| ModifiedNR.sumMatrix(fx0) > 0.5)) {
				System.out.printf("\nSolution found in: %d msec \t\t <<<<<<<<<********>>>>>>>>>",(System.currentTimeMillis() - cur));
//				for (Bus b : buses) {
//					System.out.printf(
//							"\nBus index: %d \t Bus type: %s\n" + "Bus Voltage: %.4f\n" + "Bus Active Power: %.4f\n"
//									+ "Bus Reactive Power: %.4f\n"
//									+ "-------------------------------------------------",
//							buses.indexOf(b), b.type == 0 ? "PV" : b.type == 1 ? "PQ" : "DROOP", b.voltage.getValue(),
//							b.p, b.q);
//
//				}
//				for(Bus b: buses) {
//					if(b.type==2)
//						System.out.printf("\nDroop coefficients of bus %d  \t mp: %.5f nq: %.5f\n",buses.indexOf(b),b.mp,b.nq);
//					
//				}
//				System.out.println("\nFinal Mismatches: " + fx0);
//				System.out.println("Mismatch Summation: " + ModifiedNR.sumMatrix(fx0));
//				System.out.printf("\n\nPLoss: %.5f\tPLoad: %.5f\n" + "QLoss: %.5f\tQLoad: %.5f \n", -PQLossLoad.get(0),
//						PQLossLoad.get(2), -PQLossLoad.get(1), PQLossLoad.get(3));
//				System.out.printf("\n\nTotal Active Power Demand: %.5f\n" + "Total Reactive Power Demand: %.5f \n",
//						-PQLossLoad.get(0) + PQLossLoad.get(2), -PQLossLoad.get(1) + PQLossLoad.get(3));
//				System.out.printf("\nSteady state system frequency: %.5f\n", wi);
//				System.out.printf("Steady state reference voltage value: %.5f\n", buses.get(0).voltage.getValue());
				if(ModifiedNR.sumMatrix(fx0)<min_mismatch) {
					best_mp[0] = buses.get(1).mp;
					best_mp[1] = buses.get(2).mp;
					best_mp[2] = buses.get(3).mp;
					best_nq[0] = buses.get(1).nq;
					best_nq[1] = buses.get(2).nq;
					best_nq[2] = buses.get(3).nq;
					min_mismatch = ModifiedNR.sumMatrix(fx0);
				}
				//break outerloop;
			}

			
		}
		

			System.out.printf("\nmp_1: %.5f \t mp_2: %.5f \t mp_3: %.5f",best_mp[0],best_mp[1],best_mp[2]);
			System.out.printf("\nnq_1: %.5f \t nq_2: %.5f \t nq_3 : %.5f",best_nq[0],best_nq[1],best_nq[2]);
			System.out.printf("\nMinimum mismatch (in 50 iterations): %.5f",min_mismatch);
	}
}
