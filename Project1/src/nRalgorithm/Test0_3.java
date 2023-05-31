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


public class Test0_3 {
	public static final double loadMean = 0.15;
	public static final double loadSD = loadMean * 0.1;

	public static void main(String[] args) {
		int rep = 100;
		
		NormalDistribution normal = new NormalDistribution(loadMean, loadSD);
		
		// shape 7.5 scale 3.5
		WeibullDistribution wb = new WeibullDistribution(7.5, 3.5);
		WeibullDistribution wb2 = new WeibullDistribution(7.5, 7);
		Test0_1 tst = new Test0_1();

		ModifiedNR10 mnr10 = new ModifiedNR10();
		double[][] xadmittances = mnr10.xadmittances;
		double[][] radmittances = mnr10.radmittances;

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
		boolean flag2 = false;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;

		
		int cv = 0;
		String method = "30PSO_01501_100_";
		
		

//		double[] position = {0.033231,0.023904,0.017174,0.226110,0.094650,0.085099,0.103928,0.058397}; //PSO without MGCC
		
//		double[] position = {0.0125,0.00833,0.0125,0.00833,0.1667,0.1111,0.1667,0.1111}; 
		
//		double[] position = {0.347567,	0.003632,	1.205042,	1.353550,	0.251917,	0.017686,	4.383861,	0.155639}; // PSO with MGCC
		
//		double[] position = {0.347567,	0.003632,	1.205042,	1.353550,	0.251917,	0.017686,	4.383861,	0.155639}; //PSO 30 
		
		double[] position = {-1.417865,	2.049702,	1.153322,	0.323612,	-0.722969,	0.239871,	0.598040,	0.297869,	-1.277499,	0.701896,	0.889773,	-0.045808,	
							-0.917064,	-1.990770,	-0.174837,	-1.270405,	-0.221871,	2.670507,	6.966185,	2.415106};
		double p2 = 1;
		for(int i =0; i<position.length;i++)
			position[i] = position[i] * p2 ;
			
		for (int k = 0; k < rep; k++) {
			System.out.println("\t"+k);
			int params = 60;
			int order = 2;

			X1 = null;
			fx0 = null;
			wi = 1.0;

			prev_Mismatches = Double.MAX_VALUE;
			flag = true;
			flag2 = false;
			Jacobian = null;

			boolean flag3 = true;
			while (flag3) {
				buses = new ArrayList<Bus>();
				buses.add(new Bus(1, 0, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 2, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 4, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(2, 6, params, 0, 0, position[0], position[10], 0.6, 0.8));
				buses.add(new Bus(1, 8, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 10, params, 0, 0, position[1], position[11], 0.9, 1.2));
				buses.add(new Bus(1, 12, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 14, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 16, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(2, 18, params, 0, 0, position[2], position[12], 0.6, 0.8));
				buses.add(new Bus(1, 20, params, ModifiedNR8.generateFromWind(wb2.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 22, params, 0, 0, position[3], position[13], 0.9, 1.2));
				buses.add(new Bus(1, 24, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 26, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 28, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(2, 30, params, 0, 0, position[4], position[14], 0.6, 0.8));
				buses.add(new Bus(1, 32, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 34, params, 0, 0, position[5], position[15], 0.9, 1.2));
				buses.add(new Bus(1, 36, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 38, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 40, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(2, 42, params, 0, 0, position[6], position[16], 0.6, 0.8));
				buses.add(new Bus(1, 44, params, ModifiedNR8.generateFromWind(wb2.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 46, params, 0, 0, position[7], position[17], 0.9, 1.2));
				buses.add(new Bus(1, 48, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 50, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(1, 52, params,  Math.min(-normal.sample(),0)));
				buses.add(new Bus(2, 54, params, 0, 0, position[8], position[18], 0.6, 0.8));
				buses.add(new Bus(1, 56, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 58, params, 0, 0, position[9], position[19], 0.9, 1.2));
				flag3 = ModifiedNR.checkAdequacy(buses, 12.5);
			}
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

				if (ModifiedNR.sumMatrix(fx0) < 1E-4) {

					cv++;
					flag2 = true;
					flag = false;

				} else if (prev_Mismatches + 1e1 <= ModifiedNR.sumMatrix(fx0)) {

					flag = false;
					
				} else {

					prev_Mismatches = ModifiedNR.sumMatrix(fx0);

					wi = X1.getEntry(X1.getRowDimension() - 2, 0);

					buses.get(0).voltage = new DerivativeStructure(params, order, 1,
							X1.getEntry(X1.getRowDimension() - 1, 0));

					ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
				}

			}

			if (!flag && flag2) {

				double fitness = ModifiedNR.evaluate3(buses, PQLossLoad, wi, position, fx0);

				if (ModifiedNR.sumMatrix(fx0) > 1E-4)
					fitness += ModifiedNR.sumMatrix(fx0);

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

		tst.printTo(fits, freqs, position, bV, method, bP, bQ, bPLoss, bQLoss);
		System.out.println(cv);

	}

	public void printTo(ArrayList<Double> fits, ArrayList<Double> freqs, double[] position,
			ArrayList<ArrayList<Double>> bV, String method, ArrayList<ArrayList<Double>> bP,
			ArrayList<ArrayList<Double>> bQ, ArrayList<Double> bPLoss, ArrayList<Double> bQLoss) {

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

			printWriter.print("\nReplication\tFitnessValue\tFrequency");

			for (int i = 0; i < bV.get(0).size(); i++)
				sb.append("\tBus" + (i + 1));

			for (int i = 0; i < bP.get(0).size(); i++)
				sb.append("\tBus" + (i + 1) + "_P");

			for (int i = 0; i < bQ.get(0).size(); i++)
				sb.append("\tBus" + (i + 1) + "_Q");

			printWriter.print(sb.toString());
			printWriter.print("\tPLoss\tQLoss");

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
				printWriter.printf("\t%.5f\t%.5f", bPLoss.get(fits.indexOf(d)), bQLoss.get(fits.indexOf(d)));

			}

			myWriter.close();
			System.out.println("Successfully wrote to the file.");
		} catch (IOException e) {
			System.out.println("An error occurred.");
			e.printStackTrace();
		}
	}

}
