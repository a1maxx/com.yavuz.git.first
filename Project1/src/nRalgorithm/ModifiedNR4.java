package nRalgorithm;

import java.io.File;
import java.io.FileNotFoundException;
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

public class ModifiedNR4 {

	public ArrayList<Bus> createBuses(double[] mp, double[] nq) {
		int params = 12;
	    NormalDistribution n2 = new NormalDistribution(0.001,0.0001);
		ArrayList<Bus> buses = new ArrayList<Bus>();
		buses.add(new Bus(1, 0, params, -0.0525369, -0.02188417));

		// Droop Bus 9.4E-4 0.0013
		buses.add(new Bus(2, 2, params, 0, 0, mp[0], nq[0], 1.2));
		buses.add(new Bus(2, 4, params, 0, 0, mp[1], nq[1], 1.0));
		buses.add(new Bus(2, 6, params, 0, 0, mp[2], nq[2], 1.4));

		// PQ Bus
		buses.add(new Bus(1, 8, params, n2.sample(), 0.0));
		buses.add(new Bus(1, 10, params, -0.04417614, -0.01281924));

		return buses;
	}

	public double runMNR(ArrayList<Bus> buses) {
		double wi = 1.0;
		int params = buses.size() * 2;
		int order = 2;
		int N = buses.size();

		double[][] xadmittances = new double[N][N];
		double[][] theta = new double[N][N];
		double[][] radmittances = new double[N][N];
		
		File file = new File("test2.txt");
		
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

		double[][] Jacobian = null;
		
		RealMatrix fx0 = null;
		ArrayList<Double> PQLossLoad = null;
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);
		RealMatrix X1 = null;
	
		double maxv = 0.0;
		double tolerance = 0.1;
		double sum = 0;
		boolean flag2 = false;
		double [] pd = new double[21];
		double [] qd = new double [21];
		double [] miss = new double [pd.length];
// 	    NormalDistribution n = new NormalDistribution(0.02,0.001);
		NormalDistribution n = new NormalDistribution(1.0,0.1);
		double [] PF = new double[pd.length];
		for (int f = 0; f < PF.length; f++) {
			if (f < pd.length/3)
				PF[f] = 0.4;
			else if (f < 2*pd.length/3)
				PF[f] = 0.6;
			else if (f < pd.length)
				PF[f] = 1.0;
		}
		for(int i=0 ; i<pd.length;i++) {
			pd[i]= -n.sample();
			qd[i]= pd[i]*0.1;
		}
		int p=0;
		for (int k = 0; k < pd.length/2; k++) {
			double prev_Mismatches = Double.MAX_VALUE;
			
			RealMatrix prevX1 = new Array2DRowRealMatrix(12, 1);
			RealMatrix prevfx0 = null;
			boolean flag = true;
			for(Bus b:buses) {
				if(b.type==1 && b.index!=8) {
					b.nominal_p= pd[p];
					b.nominal_q= qd[p];
					p++;
				}
			}
			
			innerloop: for (int i = 1; i <= 5 && flag; i++) {
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
					try {

						X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
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

				if (prev_Mismatches + tolerance < ModifiedNR.sumMatrix(fx0)) {

					X1 = prevX1;
					flag = false;

				} else if (ModifiedNR.sumMatrix(fx0) < 1E-5) {

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
			for (int i = 0; i < 2; i++)
				sum += Math.abs(PQLossLoad.get(i)) * 5;

			sum += findMax(buses) * 5;

			sum += Math.abs(1.0 - X1.getEntry(X1.getRowDimension() - 2, 0)) * 5 ;
			
			flag2=false;
			for (Bus b : buses) {
				if ((b.type == 2 && (b.mp < 0 || b.nq < 0 || b.nq > 1 || b.mp >= 1))) {
					if (b.nq <= 0)
						maxv = -(b.nq-1);
					if (b.nq >= 1)
						maxv = b.nq;
					if (b.mp <= 0)
						maxv = -(b.mp-1);
					if (b.mp >= 1)
						maxv = b.mp;

					flag2 = true;
				}
			}
			miss[k] = ModifiedNR.sumMatrix(fx0);
		}

		if (flag2)
			return 10000*(maxv*10e5);
		else if(findMax(miss)>10.0)
			return ModifiedNR.sumMatrix(fx0)*10000;
		else
			return sum/pd.length;

	}	
	
	static double findMax(ArrayList<Bus> buses) {
		double max = 0;
		for (Bus b : buses) {
			max = max < Math.abs(1.0 - b.voltage.getValue()) ? Math.abs(1.0 - b.voltage.getValue()) : max;

		}
		return max;
	}

	
	static double findMax(double [] arr) {
		double max = 0;
		for (Double b : arr) {
			max = max <b?b:max;

		}
		return max;
		
	}
	
	static void updateDemand(ArrayList<Bus> buses) {
		for (Bus b : buses) {
			if (b.type == 1) {
				double new_p = new NormalDistribution(0.1, 001).sample();
				b.nominal_p = new_p;
				b.p = new_p;
				double new_q = new NormalDistribution(0.05, 001).sample();
				b.nominal_q = new_q;
				b.q = new_q;
			}
		}
	}






}