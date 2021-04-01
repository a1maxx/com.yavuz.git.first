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

public class ModifiedNR3_1 {

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
		buses.add(new Bus(2, 2, theta.length * 2, 0, 0, .83216, .87445, 1.2));
		buses.add(new Bus(2, 4, theta.length * 2, 0, 0, .63259, .20593, 1.0));
		buses.add(new Bus(2, 6, theta.length * 2, 0, 0, .61869, .90717, 1.4));

		NormalDistribution n2 = new NormalDistribution(0.01, 0.001);
		// PQ Bus
		buses.add(new Bus(1, 8, params, n2.sample(), 0.0));
		buses.add(new Bus(1, 10, params, -0.05417614, -0.01281924));

	setDroops(new double[] {0.17670407, 0.83019276, 0.12948978},
				new double[] {0.84501671, 0.62637308, 0.76754906}, buses);
//		setDroops(new double[] {0.1388, 0.1388, 0.1388},
//				new double[] {  0.972, 0.972, 0.972 }, buses);

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);
		RealMatrix X1 = null;
		int order = 2;
		RealMatrix fx0 = null;
		double[][] Jacobian = null;
		double wi = 1.0;

		double[] pd = new double[720];
		double[] qd = new double[720];
		NormalDistribution n = new NormalDistribution(0.02, 0.001);
		double [] PF = new double[720];
		for (int f = 0; f < PF.length; f++) {
			if (f < 240)
				PF[f] = 0.4;
			else if (f < 480)
				PF[f] = 0.6;
			else if (f < 720)
				PF[f] = 1.0;
		}
	
		
		for(int r=1; r<20;r++) {
			StringBuilder PLosses = new StringBuilder();
			PLosses.append("\nPLosses :\t");
			StringBuilder QLosses = new StringBuilder();
			QLosses.append("\nQLosses :\t");
			StringBuilder Freqs = new StringBuilder();
			Freqs.append("\nFrequencies : \t");
			
			for (int i = 0; i < pd.length; i++) {
				pd[i] = -n.sample() * PF[i];
				qd[i] = pd[i] * 0.1;
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

			innerloop: for (int i = 1; i <= 50 && flag; i++) {
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

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {

					X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
					if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
						System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
						break innerloop;
					}
				} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

					e.printStackTrace();

				}
				if (prev_Mismatches <= ModifiedNR.sumMatrix(fx0)) {
					System.out.printf("\nMissed the local optimum. Exiting...");
					System.out.printf("\nRejected sum of mismatches: %.5f", ModifiedNR.sumMatrix(fx0));
					X1 = prevX1;
					flag = false;

				} else if (ModifiedNR.sumMatrix(fx0) < 1E-4) {
					System.out.printf("\nFull convergence achieved. Exiting...");
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

			if (!Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {

				PLosses.append("\t"+ -PQLossLoad.get(0));
				QLosses.append("\t"+ -PQLossLoad.get(1));
				Freqs.append("\t" + wi );

			}

		}
		try {
			
			myWriter.write(PLosses.toString() );
			myWriter.write(QLosses.toString());
			myWriter.write(Freqs.toString());
//			myWriter.write("\n" + "-".repeat(16) +"EOF REP "+r+"-".repeat(16));
		} catch (IOException e) {

			e.printStackTrace();
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
