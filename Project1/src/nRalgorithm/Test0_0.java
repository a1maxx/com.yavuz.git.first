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
		int rep = 100;
		
		double loadMean = 0.1;
		
		Test0 tst =  new Test0();
		
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
		ArrayList<ArrayList<Double>> bV =  new ArrayList<ArrayList<Double>>(rep);
		for(int i=0; i < rep; i++) {
		    bV.add(new ArrayList<Double>());
		}
		
		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.20);
		int cv = 0;
		String method = "CONV_";
//		double[] position = {0.374205,	0.778976,	0.007863,	0.875875,	0.343483,	0.049583}; // PSO1
//		double[] position = {0.090967,	0.620688,	0.123440,	0.432060,	0.164542,	0.45944}; // PSO2
//		double[] position = {0.212889,	0.012893,	0.145431,	0.590238,	0.377057,	0.064886}; // PSO3
//		double[] position = {0.1388,0.1388,0.1388,	0.972,	0.972,	0.972}; // oOnventional
		double[] position = {0.05551,0.05551,0.05551,0.05559,0.05559,0.05559};
		
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
			buses.add(new Bus(1, 2, params, ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
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
		


		tst.printTo(fits, freqs, position,bV,method);
		System.out.println(cv);

	}

	public void printTo(ArrayList<Double> fits, ArrayList<Double> freqs, double[] position,ArrayList<ArrayList<Double>> bV,String method) {

		  String fileName = new SimpleDateFormat("yyyyMMddHHmmss'.txt'").format(new Date());
		   try {
		    	
		    	  fileName = "zoutput_" + method + fileName  ;
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
			StringBuilder sb =  new StringBuilder();

			printWriter.print("\n Replication \t FitnessValue \t Frequency");
			
			for(int i=0;i<bV.get(0).size();i++)
				sb.append("\tBus"+(i+1));
			
			printWriter.print(sb.toString());	
		
			
			for(Double d: fits) {
				printWriter.printf("\n%s %.5f\t%.5f","Replication" + fits.indexOf(d)+"\t" ,d,freqs.get(fits.indexOf(d)));
				for(Double f:bV.get(fits.indexOf(d))) {
					printWriter.printf("\t%.5f ",f);
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
