package nRalgorithm;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

public class Test {	

	public static void main(String[] args) {

		MultivariateFunction p = new MultivariateFunction() {
			public double value(double[] x) {
				return FastMath.pow(x[0], 2) + FastMath.pow(x[1], 2);
			}

		};


		// Ybus elements
		double Y11, Y12, Y13, Y21, Y22, Y23, Y31, Y32, Y33;
		Y11 = -14;
		Y12 = 10;
		Y13 = 4;
		Y21 = 10;
		Y22 = 15;
		Y23 = 5;
		Y31 = 4;
		Y32 = 5;
		Y33 = 5;

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


		RealMatrix X0 = new Array2DRowRealMatrix();
		double[][] admittances = { { Y11, Y12, Y13 }, { Y21, Y22, Y23 }, { Y31, Y32, Y33 } };
		double[][] theta = { { t11, t12, t13 }, { t21, t22, t23 }, { t31, t32, t33 } };

		ArrayList<Bus> buses = new ArrayList<Bus>();
		buses.add(new Bus(3, admittances, theta, 0, 0, 0));
		buses.add(new Bus(0, admittances, theta, 2, -0.9, -0.5));
		buses.add(new Bus(1, admittances, theta, 4, 0.6, 0));
		buses.get(0).voltage = new DerivativeStructure(params, order, 1, 1.0);
		buses.get(0).delta = new DerivativeStructure(params, order, 0, 0);
		buses.get(2).voltage = new DerivativeStructure(params, order, 5, 1.01);

		DerivativeStructure[] equations = Bus.createEqauations2(buses);

		for (int i = 1; i <= 5; i++) {
			equations = Bus.createEqauations2(buses);
			double[] unknowns0 = {buses.get(1).delta.getValue(), buses.get(2).delta.getValue(),
					buses.get(1).voltage.getValue() };
			X0 = new Array2DRowRealMatrix(unknowns0);
			double[][] jacobs = {
					{ equations[0].getPartialDerivative(0, 0, 1, 0, 0, 0),
							equations[0].getPartialDerivative(0, 0, 0, 0, 1, 0),
							equations[0].getPartialDerivative(0, 0, 0, 1, 0, 0) },
					{ equations[2].getPartialDerivative(0, 0, 1, 0, 0, 0),
							equations[2].getPartialDerivative(0, 0, 0, 0, 1, 0),
							equations[2].getPartialDerivative(0, 0, 0, 1, 0, 0) },
					{ equations[1].getPartialDerivative(0, 0, 1, 0, 0, 0),
							equations[1].getPartialDerivative(0, 0, 0, 0, 1, 0),
							equations[1].getPartialDerivative(0, 0, 0, 1, 0, 0) } };
			RealMatrix J = new Array2DRowRealMatrix(jacobs);
			double[] functions0 = { equations[0].getValue(), equations[2].getValue(), equations[1].getValue() };
			RealMatrix fx0 = new Array2DRowRealMatrix(functions0);

			RealMatrix JI = MatrixUtils.inverse(J);
			X0 = X0.subtract(JI.multiply(fx0));
			buses.get(1).delta = new DerivativeStructure(params, order, 2, X0.getEntry(0, 0));
			buses.get(1).voltage = new DerivativeStructure(params, order, 3, X0.getEntry(2, 0));
			buses.get(2).delta = new DerivativeStructure(params, order, 4, X0.getEntry(1, 0));
			System.out.printf("\td2 = %7.6f \t iteration %d %n", X0.getEntry(0, 0), i);
			System.out.printf("\td3 = %7.6f \t iteration %d %n", X0.getEntry(1, 0), i);
			System.out.printf("\tv2 = %7.6f \t iteration %d %n", X0.getEntry(2, 0), i);
			System.out.printf("%s%n", "--------------------------------------------");
		}

	}

}
