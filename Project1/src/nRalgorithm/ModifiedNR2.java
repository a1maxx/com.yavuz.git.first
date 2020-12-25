package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class ModifiedNR2 {

	public static void main(String[] args) {

		double pi = Math.PI;

		Random random = new Random();

		double[][] xadmittances = new double[6][6];
		double[][] theta = new double[6][6];
		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		ArrayList<Double> PQLossLoad = null;

		for (int i = 0; i < xadmittances.length; i++) {
			for (int j = 0; j < xadmittances[1].length; j++) {
				xadmittances[i][j] = i == j ? 0 : (random.nextDouble() + 0.5);
				theta[i][j] = i == j ? -pi / 2 : pi / 2;
			}
		}

		for (int i = 0, len = radmittances.length; i < len; i++)
			Arrays.fill(radmittances[i], (random.nextDouble()+1)*0.2);

		ArrayList<Bus> buses = new ArrayList<Bus>();
		// PV BUS
		buses.add(new Bus(0, 0, theta.length * 2, 0, 0, 1.1));
		// Droop Bus
		buses.add(new Bus(2, 2, theta.length * 2, 0, 0, 0.0692791530265256, 0.40898603960680735));
		buses.add(new Bus(2, 4, theta.length * 2, 0, 0, -0.0692791530265256, 0.40898603960680735));
		buses.add(new Bus(2, 6, theta.length * 2, 0, 0, 0.0692791530265256, 0.40898603960680735));
		// PQ Bus
		buses.add(new Bus(1, 8, theta.length * 2, -0.8, -0.1));
		buses.add(new Bus(1, 10, theta.length * 2, -0.5, -0.2));

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

		int params = buses.size() * 2;
		int order = 2;

		double wi = 1;
		long cur = System.currentTimeMillis();
		
		for (int i = 1; i < 50; i++) {
			double w0 = 1.0;
			double v0 = 1.1;

			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix(radmittances, xadmittances, wi);

			ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

			RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi, buses.get(0).voltage.getValue());

			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations4(buses,
					Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances));

			PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

			double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, wi, w0, v0, pq, PQLossLoad, indexes);

			RealMatrix fx0 = new Array2DRowRealMatrix(mismatches);

			double[][] Jacobian = ModifiedNR.constructJacabian2(deltaVoltageOrders, pq, buses, wi, cAdmittances,
					indexes, Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
					radmittances, xadmittances);

			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

			RealMatrix X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));

			wi = X1.getEntry(X1.getRowDimension() - 2, 0);

			buses.get(0).voltage = new DerivativeStructure(params, order, 1, X1.getEntry(X1.getRowDimension() - 1, 0));

			ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, indexes, params, order);

			System.out.println("X0=\t" + X0);
			System.out.println("X1=\t" + X1);
			System.out.println("fx0=\t" + fx0);

			for (int j = 0; j < X1.getRowDimension(); j++)
				System.out.printf("\t%s = %7.6f \t iteration %d %n", "Row".concat("" + j), X1.getEntry(j, 0), i);

			System.out.printf("%s%n", "------------------------------------------------");

		}
		System.out.println("Total Time Elapsed (in msec) : " + (System.currentTimeMillis() - cur));
		for (Bus b : buses) {
			System.out.printf(
					"\nBus index: %d \t Bus type: %s\n" + "Bus Voltage= %.4f\n" + "Bus Active Power: %.4f\n"
							+ "Bus Reactive Power: %.4f\n"
							+ "-------------------------------------------------",
					buses.indexOf(b), b.type == 0 ? "PV" : b.type == 1 ? "PQ" : "DROOP", b.voltage.getValue(), b.p,
					b.q);

		}

		System.out.printf("\n\nPLoss: %.5f\tPLoad: %.5f\n"
				+ "QLoss: %.5f\tQLoad: %.5f \n", -PQLossLoad.get(0),
				PQLossLoad.get(2), -PQLossLoad.get(1), PQLossLoad.get(3));
	}

}
