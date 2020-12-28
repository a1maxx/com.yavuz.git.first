package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;

public class MicrogridFitnessFunction implements FitnessFunction{

	@Override
	public double getFitness(double[] particlePosition) {
		
		double[] mismatches = null;
		ArrayList<Double> PQLossLoad = null;
		RealMatrix X1 =null;
		
		double X11, X12, X13, X21, X22, X23, X31, X32, X33;
		X11 = 0;
		X12 = 0.1;
		X13 = 1.25;
		X21 = 0.1;
		X22 = 0;
		X23 = 0.2;
		X31 = 1.25;
		X32 = 0.2;
		X33 = 0;

		double pi = Math.PI;
		double t11, t12, t13, t21, t22, t23, t31, t32, t33;
		t11 = -pi / 2;
		t12 = pi / 2;
		t13 = pi / 2;
		t21 = pi / 2;
		t22 = -pi / 2;
		t23 = pi / 2;
		t31 = pi / 2;
		t32 = pi / 2;
		t33 = -pi / 2;

		// Known variables

		int params = 6;
		int order = 2;

		double[][] xadmittances = { { X11, X12, X13 }, { X21, X22, X23 }, { X31, X32, X33 } };
		double[][] theta = { { t11, t12, t13 }, { t21, t22, t23 }, { t31, t32, t33 } };

		double[][] radmittances = new double[xadmittances.length][xadmittances[1].length];
		for (int i = 0, len = radmittances.length; i < len; i++)
			Arrays.fill(radmittances[i], 0);

		ArrayList<Bus> buses = new ArrayList<Bus>();

		buses.add(new Bus(0, xadmittances, theta, 0, 0, 0, 0, 0));
		buses.add(new Bus(2, xadmittances, theta, 2, 0, 0, particlePosition[0], particlePosition[1]));
		buses.add(new Bus(1, xadmittances, theta, 4, -0.9, -0.8));
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

	
		double wi = 1;
		
		for (int i = 1; i < 10; i++) {
			double w0 = 1.0;
			double v0 = 1.1;
			
			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix(radmittances, xadmittances, wi);
			
			ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);
			
			RealMatrix X0  = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,buses.get(0).voltage.getValue());
			
			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations4(buses,
					Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances));
			
			PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);
			
			mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad,indexes);

			RealMatrix fx0 = new Array2DRowRealMatrix(mismatches);
			
			double[][] Jacobian = ModifiedNR.constructJacabian2(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
					Admittance.createMadmittance(cAdmittances),Admittance.createTadmittance(cAdmittances),radmittances,xadmittances);
			
			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);			
			
			
			try {
				X1 = X0.subtract(MatrixUtils.inverse(JJ).multiply(fx0));
			} catch (NonSquareMatrixException e) {
				System.err.printf("Error occured with values %d\t %d\t\n",particlePosition[0],particlePosition[1]);
				return Double.MAX_VALUE;
			} catch (MatrixDimensionMismatchException e) {
				System.err.printf("Error occured with values %f\t %f\t\n",particlePosition[0],particlePosition[1]);
				return Double.MAX_VALUE;
			} catch (DimensionMismatchException e) {
				System.err.printf("Error occured with values %f\t %f\t\n",particlePosition[0],particlePosition[1]);
				return Double.MAX_VALUE;
			} catch (NullArgumentException e) {
				System.err.printf("Error occured with values %f\t %f\t\n",particlePosition[0],particlePosition[1]);
				return Double.MAX_VALUE;
			} catch (SingularMatrixException e) {
				System.err.printf("Error occured with values %f\t %f\t\n",particlePosition[0],particlePosition[1]);
				return Double.MAX_VALUE;
			}

			wi = X1.getEntry(X1.getRowDimension() - 2, 0);

			buses.get(0).voltage = new DerivativeStructure(params, order, 1, X1.getEntry(X1.getRowDimension()-1, 0));
			
			ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, indexes, params, order);
		


		}
		double sum =0;
		for(int i=0;i<mismatches.length;i++)
			sum += mismatches[i];
		
		double loss =0;
		for(int i=0;i<PQLossLoad.size();i++)
			loss += Math.abs(PQLossLoad.get(i));
		
		sum += 4*Math.abs(1.1-X1.getEntry(X1.getRowDimension()-1, 0));
		sum += 4*Math.abs(1-X1.getEntry(X1.getRowDimension()-2, 0));
		
		return Math.abs(sum+loss);
	}

}
