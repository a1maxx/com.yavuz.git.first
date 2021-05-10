package nRalgorithm;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class Test0_0 {

	public static void main(String[] args) {
		int rep = 25;

		double loadMean = 0.25;

		Test0_0 tst = new Test0_0();

		ModifiedNR8 mnr8 = new ModifiedNR8();
		double[][] xadmittances = mnr8.xadmittances;
		double[][] radmittances = mnr8.radmittances;

		WeibullDistribution wb = new WeibullDistribution(7.5, 3.5);

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;

		ArrayList<Double> fits = new ArrayList<Double>();
		ArrayList<Double> freqs = new ArrayList<Double>();
		ArrayList<ArrayList<Double>> bV = new ArrayList<ArrayList<Double>>(rep);
		ArrayList<ArrayList<Double>> bP = new ArrayList<ArrayList<Double>>(rep);
		ArrayList<ArrayList<Double>> bQ = new ArrayList<ArrayList<Double>>(rep);
		ArrayList<Double> bPLoss = new ArrayList<Double>();
		ArrayList<Double> bQLoss = new ArrayList<Double>();

		for (int i = 0; i < rep; i++) {
			bV.add(new ArrayList<Double>());
			bP.add(new ArrayList<Double>());
			bQ.add(new ArrayList<Double>());
		}

		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.20);
		int cv = 0;
		String method = "PSO1_";
//		double[] position = { 0.000157, 0.439875, 0.828400, 0.017292 }; // PSO1
		double[] position = {0.001141,	0.705643,	0.085975,	0.380136}; // PSO1__iter80_npar10_budget50
//		double[] position = {0.659190,	0.000189,	0.003205,	0.108217};// PSO2__iter100_npar30_budget50
		

//		double[] position = {0.004145,	0.008543,	0.056562,	0.458590}; // PSO2__iter80_npar10_budget50
//		double[] position = {0.000185,	0.178681,	0.005492,	0.366317};// PSO2__iter100_npar30_budget50
//		double[] position = {0.003576,	0.618300,	0.013196,	0.102596};// PSO2_npar10_niter100_n10_budget150
		
		
//		double[] position = {0.000594,	0.414312,	0.389496,	0.005235}; // PSO3__iter80_npar10_budget50
//		double[] position = {0.761462,	0.000132,	0.722829,	0.001252}; // PSO3__iter100_npar30_budget50
//		double[] position = {0.896687,	0.000196,	0.352531,	0.011333}; // PSO3_npar10_niter100_n10_budget150


//		double[] position = {0.1679167,0.1119444,0.1654167,0.1102778}; // Conventional_newSetup

		for (int k = 0; k < rep; k++) {
			int params = 12;
			int order = 2;

			X1 = null;
			fx0 = null;
			wi = 1.0;

			prev_Mismatches = Double.MAX_VALUE;
			flag = true;

			Jacobian = null;

			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, params, -normal.sample()));
			buses.add(new Bus(1, 2, params, -normal.sample()));
			buses.add(new Bus(1, 4, params, -normal.sample()));
			buses.add(new Bus(2, 6, params, 0, 0, position[0], position[2], 0.8));
			buses.add(new Bus(1, 8, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
			buses.add(new Bus(2, 10, params, 0, 0, position[1], position[3], 1.2));

			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 1000 && flag; i++) {
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

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

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

				if (ModifiedNR.sumMatrix(fx0) < 1E-3) {

					cv++;

					flag = false;

				} else if (prev_Mismatches + 1e2 <= ModifiedNR.sumMatrix(fx0)) {

					flag = false;
				} else {

					prev_Mismatches = ModifiedNR.sumMatrix(fx0);

					wi = X1.getEntry(X1.getRowDimension() - 2, 0);

					buses.get(0).voltage = new DerivativeStructure(params, order, 1,
							X1.getEntry(X1.getRowDimension() - 1, 0));

					ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
				}

			}

			if (!flag) {
				double fitness = Math.abs(PQLossLoad.get(0)) + Math.abs(PQLossLoad.get(1));
				for (Bus b : buses)
					fitness = fitness + Math.abs(1.0 - b.voltage.getValue());

				if (ModifiedNR.sumMatrix(fx0) > 1E-4)
					fitness += ModifiedNR.sumMatrix(fx0);

				fitness += Math.abs(1 - wi);
				fits.add(fitness);
				freqs.add(wi);
				bPLoss.add(Math.abs(PQLossLoad.get(0)));
				bQLoss.add(Math.abs(PQLossLoad.get(1)));

				for (Bus b : buses) {
					bV.get(k).add(b.voltage.getValue());
					bP.get(k).add(b.p);
					bQ.get(k).add(b.q);
				}

			}

		}

		tst.printTo(fits, freqs, position, bV, method,bP,bQ);
		System.out.println(cv);

	}

	public void printTo(ArrayList<Double> fits, ArrayList<Double> freqs, double[] position,
			ArrayList<ArrayList<Double>> bV, String method, ArrayList<ArrayList<Double>> bP,
			ArrayList<ArrayList<Double>> bQ) {

		String fileName = new SimpleDateFormat("yyyyMMddHHmmss'.txt'").format(new Date());
		try {

			fileName = "zoutput_" + method + fileName;
			File myObj = new File(fileName);
			if (myObj.createNewFile()) {
				System.out.println("File created: " + myObj.getName());
			} else {
				System.out.println("File already exists.");
			}
		} catch (IOException e) {
			System.out.println("An error occurred.");
			e.printStackTrace();
		}

		try {
			FileWriter myWriter = new FileWriter(fileName);
			PrintWriter printWriter = new PrintWriter(myWriter);
			StringBuilder sb = new StringBuilder();

			printWriter.print("\n Replication \t FitnessValue \t Frequency");

			for (int i = 0; i < bV.get(0).size(); i++)
				sb.append("\tBus" + (i + 1));
			
			for (int i = 0; i < bP.get(0).size(); i++)
				sb.append("\tBus" + (i + 1)+"_P");
			
			for (int i = 0; i < bQ.get(0).size(); i++)
				sb.append("\tBus" + (i + 1)+"_Q");
			
			printWriter.print(sb.toString());

			for (Double d : fits) {
				printWriter.printf("\n%s %.5f\t%.5f", "Replication" + fits.indexOf(d) + "\t", d,
						freqs.get(fits.indexOf(d)));
				for (Double f : bV.get(fits.indexOf(d))) {
					printWriter.printf("\t%.5f ", f);
				}
				for (Double f : bP.get(fits.indexOf(d))) {
					printWriter.printf("\t%.5f ", f);
				}
				for (Double f : bQ.get(fits.indexOf(d))) {
					printWriter.printf("\t%.5f ", f);
				}
			}

			myWriter.close();
			System.out.println("Successfully wrote to the file.");
		} catch (IOException e) {
			System.out.println("An error occurred.");
			e.printStackTrace();
		}
	}

}
