package nRalgorithm;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class Test0 {

	public static void main(String[] args) {
		int rep = 10;
		double loadMean = 0.05;
		Test0 tst =  new Test0();
		
		ModifiedNR8 mnr8 = new ModifiedNR8();
		double[][] xadmittances = mnr8.xadmittances;
		double[][] radmittances = mnr8.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;

		ArrayList<Double> fits = new ArrayList<Double>();
		ArrayList<Double> freqs = new ArrayList<Double>();
		ArrayList<ArrayList<Double>> bV =  new ArrayList<ArrayList<Double>>(rep);
		for(int i=0; i < rep; i++) {
		    bV.add(new ArrayList<Double>());
		}
		
		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.1);
		int cv = 0;
//		double[] position = { 0.019701, 0.474667, 0.793566, 0.030766, 0.535109, 0.476926 };
//		double[] position = { 0.610640,	1.741969,	-1.643483, 0.904923,	1.968637,	2.187067};
//		double[] position = {1.055217,	-0.575010,	0.028957,	1.351441,	0.602383,	0.956652};
		double[] position = {-0.1388,-0.1388,-0.1388,	0.972,	0.972,	0.972};
		
		for (int k = 0; k < rep; k++) {
			int params = 12;
			int order = 2;
			
			X1 = null;
			fx0 = null;
			wi = 1.0;

			prev_Mismatches = Double.MAX_VALUE;
			flag = true;

			Jacobian = null;
			double PF = 0.328;

			

			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, params, -normal.sample()));
			buses.add(new Bus(1, 2, params, 0.05, 0.05 * PF));
			buses.add(new Bus(1, 4, params, -normal.sample()));
			buses.add(new Bus(2, 6, params, 0, 0, position[0], position[3], 3.0));
			buses.add(new Bus(2, 8, params, 0, 0, position[1], position[4], 3.0));
			buses.add(new Bus(2, 10, params, 0, 0, position[2], position[5], 3.0));

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
			
				
				for(Bus b:buses) {
					bV.get(k).add(b.voltage.getValue());
					
				}
				
				
			}

		}
		


		tst.printTo(fits, freqs, position,bV);
		System.out.println(cv);

	}

	public void printTo(ArrayList<Double> fits, ArrayList<Double> freqs, double[] position,ArrayList<ArrayList<Double>> bV ) {

		try {
			FileWriter myWriter = new FileWriter("output.txt");
			PrintWriter printWriter = new PrintWriter(myWriter);

			printWriter.print("Droop Coefficients mp and nq:\t");
			for (Double d : position) 
				printWriter.printf("%.5f\t", d);
			printWriter.print("\nFitness Value \t Frequency");
			
			for(Double d: fits) {
				printWriter.printf("\n%.5f\t%.5f\n",d,freqs.get(fits.indexOf(d)));
				for(Double f:bV.get(fits.indexOf(d))) {
					printWriter.printf("\nBus %d voltage: \t%.5f ",bV.get(fits.indexOf(d)).indexOf(f),f);
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
