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

public class ModifiedNR0 {

	public static void main(String[] args) {
		
		int nofB = 6;
		double[][] xadmittances = new double[nofB][nofB];
		double[][] theta = new double[nofB][nofB];
	
		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		ArrayList<Double> PQLossLoad = null;
		ArrayList<ArrayList<DerivativeStructure>> pq= null;
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
		boolean flag2 = true;
		double wi = 1.0;
		RealMatrix prevX1 = new Array2DRowRealMatrix(6, 1);
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;
		Random random = new Random();
		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		long cur = System.currentTimeMillis();
		double[][] Jacobian = null;

		int params = 12;
		int order = 2;
		X1 = null;
		fx0 = null;
		wi = 1.0;
		prevfx0 = null;
		prevX1 = new Array2DRowRealMatrix(params, 1);
		prev_Mismatches = Double.MAX_VALUE;
		flag = true;

		Jacobian = null;
		
			flag = true;
			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, params, -0.04832, -0.07588042));
			buses.add(new Bus(1, 2, params, -0.0, -0.0));
			buses.add(new Bus(1, 4, params, -0.057376019, -0.093323308));
//			buses.add(new Bus(2, 6, params, 0, 0, -0.014902, 0.613651, 3.0));
//			buses.add(new Bus(2, 8, params, 0, 0, -0.019183, 0.641440, 3.0));
//			buses.add(new Bus(2, 10, params, 0, 0, -0.007715, 0.550888, 3.0));
			buses.add(new Bus(2, 6, params, 0, 0, -0.014902, 0.613651, 3.0));
			buses.add(new Bus(2, 8, params, 0, 0, -0.019183, 0.641440, 3.0));
			buses.add(new Bus(2, 10, params, 0, 0, -0.007715, 0.550888, 3.0));
			ModifiedNR3.setDroops(new double[] {0.610640,	1.741969,	-1.643483},
					new double[] {  0.904923,	1.968637,	2.187067}, buses);
			
			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 1000 && flag; i++) {
				double w0 = 1.00;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

				PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);
 
				RealMatrix X0 = ModifiedNR.setUnknownValues2(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

				RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

				RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

				pq = ModifiedNR.createEquations5(buses, mAdmittances,
						tAdmittances);

				double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

				fx0 = new Array2DRowRealMatrix(mismatches);

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {
					try {
						
//						X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
						
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

					e.printStackTrace();

					break innerloop;

				}
				if (ModifiedNR.sumMatrix(fx0) < 1E-8) {
					boolean flag3 = true;
				
					if (flag3) {
						System.out.println("\nFull convergence achieved. Exiting...");
						flag = false;
						flag2 = false;
					}

				} else if (prev_Mismatches + 1e2 <= ModifiedNR.sumMatrix(fx0)) {

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
					
					ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);
				}

			}
		
		if (!Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
			System.out.println("\nTotal Time Elapsed (in msec) : " + (System.currentTimeMillis() - cur));
			for (Bus b : buses) {
				System.out.printf(
						"\nBus index: %d \t Bus type: %s\n" + "Bus Voltage: %.5f\n" + "Bus Angle: %.5f\n"
								+ "Bus Active Power: %.5f\n" + "Bus Reactive Power: %.5f\n"
								+ "-------------------------------------------------\n",
						buses.indexOf(b), b.type == 0 ? "PV" : b.type == 1 ? "PQ" : "DROOP", b.voltage.getValue(),
						b.delta.getValue(), b.p, b.q );

			}
			
			for(Bus b: buses) {
				if(buses.indexOf(b)!=0) {
					System.out.printf("\n\nBus %d \nCalculated P: %.5f\n"
							+ "Calculated Q: %.5f", buses.indexOf(b),
							pq.get(0).get(buses.indexOf(b)-1).getValue(),
							pq.get(1).get(buses.indexOf(b)-1).getValue());
				}
				
			}

			for (Bus b : buses)
				if (b.type == 2)
					System.out.printf("\n\nDroop coefficients of bus %d  \t mp: %.5f nq: %.5f\n", buses.indexOf(b), b.mp,
							b.nq);

			// System.out.println("\nFinal Mismatches = " + prevfx0);
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
