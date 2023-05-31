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


public class Test0_1 {
	public static final double loadMean = 0.015;
	public static final double loadSD = loadMean * 0.5;

	public static void main(String[] args) {
		int rep = 10;
		
		NormalDistribution normal = new NormalDistribution(loadMean, loadSD);
		
		// shape 7.5 scale 3.5
		WeibullDistribution wb = new WeibullDistribution(7.5, 3.5);
		WeibullDistribution wb2 = new WeibullDistribution(7.5, 3.5);
		Test0_1 tst = new Test0_1();

		ModifiedNR9_1 mnr9 = new ModifiedNR9_1();
		double[][] xadmittances = mnr9.xadmittances;
		double[][] radmittances = mnr9.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.00;

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
		String method = "12PSO_01501_1000_";
		
		

//		double[] position = {0.033231,0.023904,0.017174,0.226110,0.094650,0.085099,0.103928,0.058397}; //PSO without MGCC
		
//		double[] position = {0.0125,0.0083,0.0125,0.0083,
//							 0.1667,0.1111,0.1667,0.1111};  //Honest Conventional
		
//		double ff = 0.3;
//		double[] position = {0.08958335+ff,0.0597222+ff,0.08958335+ff,0.0597222+ff,
//				0.02083335+ff,0.05138335+ff,0.02083335+ff,0.05138335+ff};  //Honest Conventional 2
		
//		double[] position = {0.347567,	0.003632,	1.205042,	1.353550,	
//							 0.251917,	0.017686,	4.383861,	0.155639}; // PSO with MGCC
		
//		double[] position = {0.152746,	-0.082997,	0.810816,	-0.001912,	
//							0.775256,	1.088017,	0.097911,	0.930320}; //
			
//		double[] position = {0.2873,  0.2892,  0.4535,  0.3096,
//				-0.7215, -0.7207, -0.9999, -0.8626}; // Scenario-based 3 scenarios
		
//		double[] position = {0.2895,  0.2942,  0.4253,  0.2849, -0.7728, -0.7808, -1.    ,
//			       -0.843 }; // Scenario-based 10 scenarios
		
//		double[] position = {0.2694,  0.2866,  0.4394,  0.2875, -0.6138, -0.6473, -0.9993,
//			       -0.6682}; // Scenario-based 20 scenarios
		
//		double[] position = {0.3046, 0.3192, 0.4516, 0.3008, 0.5114, 0.5414, 0.752 , 0.5008}; // Scenario-based 40 scenarios
//		double[] position = {0.29  , 0.2942, 0.4324, 0.2884, 0.3056, 0.3153, 0.4449, 0.3001}; 
		
		
//		double[] position = {0.3064, 0.3168, 0.4514, 0.3006, 0.285 , 0.3003, 0.4105, 0.2751}; // Scenario-based 70 scenarios
		
//		double[] position = {.2805, 0.2876, 0.445 , 0.2939, 0.2751, 0.2848, 0.4285, 0.2819}; 
		
		
//		double[] position = {0.377321,	0.205182,	0.211784,	0.326381,	
//							0.951991,	0.862956,	0.931548,	0.992576}; 
		
		
//		double[] position = {0.3002, 0.3159, 0.462 , 0.3005, 0.3   , 0.3137, 0.4296, 0.3}; /// Scenario 20 0.1SD
//		double[] position = {0.2423, 0.2574, 0.417 , 0.2712, 0.3   , 0.3009, 0.4724, 0.3084}; /// Scenario 20 0.5SD
//		
		double[] position = {0.253066, 0.277613, 0.354755, 0.232186, 0.325022, 0.35594 , 0.4076  , 0.325003}; /// Scenario 20 0.5SD

		
		
//		double [] position = {0.649955,	0.466609,	-0.014265,	-0.003010,	0.741878,	0.760559,	0.679216,	-0.040573}; ///PSO
		
//		double[] position = {0.033231,0.023904,0.017174,0.226110,0.094650,0.085099,0.103928,0.058397}; //PSO without MGCC
		
//		double[] position = {0.009905,	0.017098,0.009905,	0.017098,	0.667299,	0.545103,0.667299,	0.545103};//PSO without MGCC
		
//		double[] position = {0.377321,	0.205182,	0.211784,	0.326381,	0.951991,	0.862956,	0.931548,	0.992576};//PSO with MGCC
		

		
		for (int k = 0; k < rep; k++) {
			System.out.println("\t"+k);
			int params = 24;
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
				buses.add(new Bus(1, 0, params, Math.min(0,-normal.sample())));//1
				buses.add(new Bus(1, 2, params, Math.min(0,-normal.sample())));//2
				buses.add(new Bus(1, 4, params, Math.min(0,-normal.sample())));//3
				buses.add(new Bus(2, 6, params, 0, 0, position[0], position[4], 0.6, 0.8));//4
				buses.add(new Bus(1, 8, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));//5
				buses.add(new Bus(2, 10, params, 0, 0, position[1], position[5], 0.9, 1.2));//6
				buses.add(new Bus(1, 12, params, Math.min(0,-normal.sample())));//7
				buses.add(new Bus(1, 14, params, Math.min(0,-normal.sample())));//8
				buses.add(new Bus(2, 16, params, 0, 0, position[2], position[6], 0.6, 0.8));//9
				buses.add(new Bus(1, 18, params, Math.min(0,-normal.sample())));//10
				buses.add(new Bus(1, 20, params, ModifiedNR8.generateFromWind(wb2.sample(), 3.5, 20.0, 14.5, 0.75)));//11
				buses.add(new Bus(2, 22, params, 0, 0, position[3], position[7], 0.9, 1.2));//12

				flag3 = ModifiedNR.checkAdequacy(buses, 5);
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
//					JA-MNR
					
					RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
					double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
							radmittances, xadmittances, DELTA, X0);
					RealMatrix aj = new Array2DRowRealMatrix(arJacob);
					DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
					X1 = X0.add(DELTA);
					
//					Conventional
					
//					X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));

					if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
						System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
						break innerloop;
					}
				} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

					e.printStackTrace();
				}

				if (ModifiedNR.sumMatrix(fx0) < 1E-3) {

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
