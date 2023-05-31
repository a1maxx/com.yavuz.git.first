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
import org.apache.commons.math3.distribution.BetaDistribution;

public class Test0_5 {
	public static final double loadMean = 0.0015;
	public static final double loadSD = loadMean * 0.1;
	public static final double v0 = 1.01;
	public static final double w0 = 1.00;
	public static final double SENSITIVITY = 1E-3;
	
	public static void main(String[] args) {
	
		BetaDistribution beta = new BetaDistribution(0.4,8.56);
		NormalDistribution normal = new NormalDistribution(loadMean, loadSD);

		// shape 7.5 scale 3.5
		WeibullDistribution wb = new WeibullDistribution(7.5, 3.5);
		WeibullDistribution wb2 = new WeibullDistribution(7.5, 3.5);
		Test0_5 tst = new Test0_5();

		ModifiedNR11 mnr10 = new ModifiedNR11();
		double[][] xadmittances = mnr10.xadmittances;
		double[][] radmittances = mnr10.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;
		int rep = 1000;
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
		String method = "33PSO_001501_5002_";

//		double[] position = {1.761051,	0.177863,	1.089879,	0.888169,	0.111381,	0.916111,	
//		0.752919,	0.752796,	2.315153,	1.321335,	2.016096,	0.420088};

//		double[] position = {1.153545,	0.056984,	1.379426,	1.062268,	0.145582,	3.601845,	0.856098,	
//		1.170407,	0.847968,	0.698722,	2.162787,	0.835679}; //New Solution Long Run PSO_ 0.015,0.1

//		
//		double[] position = {-0.404628,	0.023286,	-0.029881,	0.060516,	
//								 -2.132032,	-2.340725,	-3.909333,	-3.104619};  // 11
////		
//		double[] position = {-0.764751,	-0.041600,	0.087224,	-0.044049,	
//								-1.377096,	-3.257998,	-5.525477,	-1.790784}; //PSO33 0.0015 0.1 250 // 22

//		double [] position = {0.8,1.9,4.5,2.8,-12.07,-19.97,-7.57,-17.94};
//		
//		double [] position = {0.8     , 1.9     , 4.5     , 2.8     , 5.82068 , 9.720742,
//	     3.814654, 9.762857};
//		
//		double [] position = {0.8     , 1.9     , 4.5     , 2.8     , 1.996364, 3.225866,
//			       2.208224, 2.012584};
//		

//		double [] position = {-0.026,	1.685,	-3.189,	-0.077,	1.836,	-1.465,	-1.504,	-1.787}; //ABC
//		
		
		
		
//		double [] position = {-0.740,	-0.439,	-1.011,	-0.274,	-1.053,	-0.657,	-5.899,	2.666}; // ABC zoutput_33PSO_001501_5002_20221113214644.txt

		double [] position = {1.943672,	1.116734,	0.033168,	-2.988463,	-10.624786,	5.916720,	-0.578030,	-7.299559}; // PSO_33_Fin 
		
		//		IPOPT WORKIN -- KEEP THIS /////////////////////////////////////////////////////////
//		double[] position = { 0.011827, 0.00870659, 0.00688067, 4.95530569, 0.00907329, ///
//				0.00820676, 1.13718583, 0.97323469, 1.14765667, 1.50637516, ///
//				0.96753198, 1.26801567 }; ///
		///////////////////////////////////////////////////////////////////////////////////
		
		for (int k = 0; k < rep; k++) {
			System.out.println("\t" + k);
			int params = 66;
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
				buses.add(new Bus(1, 0, params,  Math.min(-normal.sample(),0)));//1
				buses.add(new Bus(1, 2, params,  Math.min(-normal.sample(),0)));//2
				buses.add(new Bus(1, 4, params,  Math.min(-normal.sample(),0)));//3
				buses.add(new Bus(1, 6, params,  Math.min(-normal.sample(),0)));//4
				buses.add(new Bus(1, 8, params,  Math.min(-normal.sample(),0)));//5
				buses.add(new Bus(1, 10, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75))); // 6
				buses.add(new Bus(1, 12, params,  Math.min(-normal.sample(),0)));//7
				buses.add(new Bus(1, 14, params,  Math.min(-normal.sample(),0)));//8
				buses.add(new Bus(1, 16, params,  Math.min(-normal.sample(),0)));//9
				buses.add(new Bus(1, 18, params, ModifiedNR8.generateFromWind(wb2.sample(), 3.5, 20.0, 14.5, 0.75))); //10
				buses.add(new Bus(1, 20, params,  Math.min(-normal.sample(),0)));//11
				buses.add(new Bus(1, 22, params,  Math.min(-normal.sample(),0)));//12
				buses.add(new Bus(2, 24, params, 0, 0, position[0], position[4], 1.8, 2.4)); //13
				buses.add(new Bus(1, 26, params,  Math.min(-normal.sample(),0)));//14
				buses.add(new Bus(2, 28, params, 0, 0, position[1], position[5], 1.2, 1.6)); //15
				buses.add(new Bus(1, 30, params,  Math.min(-normal.sample(),0)));//16
				buses.add(new Bus(1, 32, params,  Math.min(-normal.sample(),0)));//17
				buses.add(new Bus(1, 34, params,  Math.min(-normal.sample(),0)));//18
				buses.add(new Bus(1, 36, params,  Math.min(-normal.sample(),0)));//19
				buses.add(new Bus(1, 38, params,  Math.min(0.02,beta.sample()*66*0.00014)));//20
				buses.add(new Bus(1, 40, params,  Math.min(-normal.sample(),0)));//21
				buses.add(new Bus(1, 42, params,  Math.min(-normal.sample(),0)));//22
				buses.add(new Bus(1, 44, params,  Math.min(-normal.sample(),0)));//23
				buses.add(new Bus(1, 46, params,  Math.min(-normal.sample(),0)));//24
				buses.add(new Bus(2, 48, params, 0, 0, position[2], position[6], 1.2, 1.6)); //25
				buses.add(new Bus(1, 50, params,  Math.min(-normal.sample(),0)));//26
				buses.add(new Bus(1, 52, params,  Math.min(-normal.sample(),0)));//27
				buses.add(new Bus(1, 54, params,  Math.min(-normal.sample(),0)));//28
				buses.add(new Bus(1, 56, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75))); //29
				buses.add(new Bus(1, 58, params,  Math.min(-normal.sample(),0)));//30
				buses.add(new Bus(1, 60, params,  Math.min(-normal.sample(),0)));//31
				buses.add(new Bus(1, 62, params,  Math.min(-normal.sample(),0)));//32
				buses.add(new Bus(2, 64, params, 0, 0, position[3], position[7], 1.2, 1.6)); //33
				

				flag3 = ModifiedNR.checkAdequacy(buses, 7);
			}
			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 1000 && flag; i++) {

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

				if (ModifiedNR.sumMatrix(fx0) < SENSITIVITY) {

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

				if (ModifiedNR.sumMatrix(fx0) > SENSITIVITY)
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
