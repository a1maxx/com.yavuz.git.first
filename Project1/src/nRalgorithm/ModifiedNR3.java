package nRalgorithm;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;

public class ModifiedNR3 {

	public static void main(String[] args) {

		double[][] xadmittances = new double[6][6];
		double[][] theta = new double[6][6];
		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		ArrayList<Double> PQLossLoad = null;
		FileWriter myWriter = null;
		try {
			myWriter = new FileWriter("output.txt");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		File file = new File("test.txt");
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
		int params = 12;

		ArrayList<Bus> buses = new ArrayList<Bus>();
		// Reference Bus (PQ)
		buses.add(new Bus(1, 0, params, -0.0525369, -0.02188417));
		// Droop Buses
		buses.add(new Bus(2, 2, theta.length * 2, 0, 0, .83216, .87445, 1.5));
		buses.add(new Bus(2, 4, theta.length * 2, 0, 0, .63259, .20593, 1.5));
		buses.add(new Bus(2, 6, theta.length * 2, 0, 0, .61869, .90717, 1.5));
		
		NormalDistribution n2 = new NormalDistribution(0.01,0.01);
		// PQ Bus
		double gen = n2.sample();
		buses.add(new Bus(1, 8, params, 0, 0));
		buses.add(new Bus(1, 10, params, 0.05417614, 0.01281924));
//
//		setDroops(new double[] {0.07771463, 0.13652186, 0.13962579},
//				new double[] {  0.98773779, 0.95290375, 0.97727852}, buses);
		setDroops(new double[] {0.1388, 0.1388, 0.1388},
				new double[] {  0.972, 0.972, 0.972 }, buses);
//		setDroops(new double[] {0.0007,0.0007, 0.0007},
//		new double[] {  0.021, 0.021, 0.021 }, buses);

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);
		RealMatrix X1 = null;
		int order = 2;
		RealMatrix fx0 = null;
		double[][] Jacobian = null;
	
		double wi = 1.0;
		long cur = System.currentTimeMillis();
		
		double[] pd = new double[6];
		double[] qd = new double[6];
		 NormalDistribution n = new NormalDistribution(0.45,0.001);

		for (int i = 0; i < pd.length; i++) {
			pd[i] = n.sample();
			qd[i] = pd[i] *0.5;
		}
		
		int p = 0;
		for (int k = 0; k < pd.length / 2; k++) {
			double prev_Mismatches = Double.MAX_VALUE;
			RealMatrix prevX1 = new Array2DRowRealMatrix(12, 1);
			RealMatrix prevfx0 = null;
			boolean flag = true;
			for (Bus b : buses) {
				if (b.type == 1 && b.index != 8) {
					b.nominal_p = pd[p];
					b.nominal_q = qd[p];
					p++;
				}
			}

			innerloop: for (int i = 1; i <= 200 && flag; i++) {
				double w0 = 1.00;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0 ,v0);

				RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

				PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

				ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances));

				double[] mismatches = ModifiedNR.calculateMismatchMatrix3(buses, pq, PQLossLoad, indexes);

				fx0 = new Array2DRowRealMatrix(mismatches);

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {
					try {

						X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
						
						//System.out.println("-".repeat(10)+ ModifiedNR.convergenceCheck(X1, X0));
						
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
				if (ModifiedNR.sumMatrix(fx0)/(X1.getRowDimension()-5) < 1E-4) {
					System.out.println("Full convergence achieved. Exiting...");
					flag = false;
	
				} else if (prev_Mismatches <= ModifiedNR.sumMatrix(fx0)) {
					for (int f = 0; f < Jacobian.length; f++) {
						for (int j = 0; j < Jacobian[1].length; j++) {
							System.out.printf("%2.2f \t", Jacobian[f][j]);
						}
						System.out.println();
					}
					
					System.out.printf("\nMissed the local optimum. Exiting... \t At iteration : %d",i);
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

//				System.out.printf("\nLatest sum of mismatches = %.5f \n\n", prev_Mismatches);
//				for (int j = 0; j < X1.getRowDimension(); j++) {
//					if (j < deltaVoltageOrders.get(0).size())
//						System.out.printf("\t%s = %7.6f \t iteration %d %n", "Row".concat("" + j),
//								X1.getEntry(j, 0) * 180 / Math.PI, i);
//					else
//						System.out.printf("\t%s = %7.6f \t iteration %d %n", "Row".concat("" + j), X1.getEntry(j, 0),
//								i);
//				}
//				System.out.printf("%s%n", "------------------------------------------------");

			}

			if (!Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
				System.out.println("\nTotal Time Elapsed (in msec) : " + (System.currentTimeMillis() - cur));
				for (Bus b : buses) {
					System.out.printf(
							"\nBus index: %d \t Bus type: %s\n" + "Bus Voltage: %.4f\n" + "Bus Angle: %.4f\n"
									+ "Bus Active Power: %.4f\n" + "Bus Reactive Power: %.4f\n"
									+ "-------------------------------------------------",
									buses.indexOf(b), b.type == 0 ? "PV" : b.type == 1 ? "PQ" : "DROOP", b.voltage.getValue(),
											b.delta.getValue() * 180.0 / Math.PI, b.p, b.q);
					
				}

				for (Bus b : buses)
					if (b.type == 2)
						System.out.printf("\nDroop coefficients of bus %d  \t mp: %.5f nq: %.5f\n", buses.indexOf(b),
								b.mp, b.nq);

				System.out.println("\nFinal Mismatches = " + prevfx0);
				System.out.println("\nMismatch Summation: " + prev_Mismatches);
				System.out.printf("\nPLoss: %.5f\tPLoad: %.5f\n" + "QLoss: %.5f\tQLoad: %.5f \n", -PQLossLoad.get(0),
						PQLossLoad.get(2), -PQLossLoad.get(1), PQLossLoad.get(3));
				System.out.printf("\nTotal Active Power Demand: %.5f\n" + "Total Reactive Power Demand: %.5f \n",
						-PQLossLoad.get(0) + PQLossLoad.get(2), -PQLossLoad.get(1) + PQLossLoad.get(3));
				try {
					myWriter.write(" " + wi);
				} catch (IOException e) {
			
					e.printStackTrace();
				}
				
	
				System.out.printf("\nSteady state system frequency: %.5f\n", wi);
				System.out.printf("Steady state reference voltage value: %.5f\n", buses.get(0).voltage.getValue());
			}

		}
		try {
			myWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	static void setDroops(double[] mp, double[] nq, ArrayList<Bus> buses) {
		int i = 0;
		for (Bus b : buses) {
			if (b.type == 2) {
				b.mp = mp[i];
				b.nq = nq[i];
				i++;
			}

		}

	}
}
