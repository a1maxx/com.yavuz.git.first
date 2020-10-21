package nRalgorithm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.analysis.*;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
//import org.apache.commons.math3.complex.Complex;

import java.lang.Math;

public class Main {

	public static void main(String[] args) {

		
		DerivativeStructure x = new DerivativeStructure(2,2,0,1);
		DerivativeStructure y = new DerivativeStructure(2,2,1,2);
		DerivativeStructure ff = x.pow(2).add(y.pow(2));
		System.out.println(ff.getValue());
		x = new DerivativeStructure(2,2,0,3);
		ff = x.pow(2).add(y.pow(2));
		System.out.println(ff.getValue());
		MultivariateFunction p =  new MultivariateFunction(){
		public double value(double [] x) {
			return FastMath.pow(x[0], 2) + FastMath.pow(x[1], 2);
		}
		
		
		};
	
		System.out.println(p.value(new double [] { 1, 2 }));
		
		// Power injections
		double p2, q2, p3;
		p2 = -0.9;
		q2 = -0.5;
		p3 = 0.6;

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
		double d2RealValue = 0.0;
		double d3RealValue = 0.0;
		double v2RealValue = 1.0;

		RealMatrix X = new Array2DRowRealMatrix();
		RealMatrix X0 = new Array2DRowRealMatrix();
		double[][] admittances = { { Y11, Y12, Y13 }, { Y21, Y22, Y23 }, { Y31, Y32, Y33 } };
		double[][] theta = { { t11, t12, t13 }, { t21, t22, t23 }, { t31, t32, t33 } };
		// RealMatrix add = new Array2DRowRealMatrix(admittances);

		Bus[] buses = new Bus[3];
		buses[0] = new Bus(3, admittances, theta, 0, 0, 0);
		buses[1] = new Bus(0, admittances, theta, 2, -0.9, -0.5);
		buses[2] = new Bus(1, admittances, theta, 4, 0.6, 0);
		buses[0].voltage = new DerivativeStructure(params, order, 1, 1.0);
		buses[0].delta = new DerivativeStructure(params, order, 0, 0);
		buses[2].voltage = new DerivativeStructure(params, order, 5, 1.01);

		DerivativeStructure d2 = new DerivativeStructure(params, order, 0, d2RealValue);
		DerivativeStructure d3 = new DerivativeStructure(params, order, 1, d3RealValue);
		DerivativeStructure v2 = new DerivativeStructure(params, order, 2, v2RealValue);
		DerivativeStructure v1 = new DerivativeStructure(params, order, 3, 1.0);
		DerivativeStructure d1 = new DerivativeStructure(params, order, 3, 0);
		DerivativeStructure v3 = new DerivativeStructure(params, order, 3, 1.01);

		DerivativeStructure fp2 = v2.multiply(v1.multiply(Y21)).multiply(d2.negate().add(d1.add(t21)).cos())
				.add(v2.pow(2).multiply(Y22 * FastMath.cos(t22)))
				.add(v2.multiply(v3.multiply(Y23)).multiply((d2.negate().add(d3.add(t23))).cos())).subtract(p2);
		DerivativeStructure fp3 = (d3.negate().add(d1.add(t31)).cos()).multiply(v1.multiply(v3.multiply(Y31)))
				.add((d3.negate().add(d2).add(t32).cos()).multiply(v2.multiply(v3.multiply(Y32))))
				.add(v3.pow(2).multiply(Y33).multiply(Math.cos(t33))).subtract(p3);
		DerivativeStructure fq2 = v2.multiply(v1.multiply(-Y21)).multiply((d2.negate().add(d1.add(t21))).sin())
				.add(v2.pow(2).multiply(-Y22).multiply(FastMath.sin(t22)))
				.add(v3.multiply(v2.multiply(-Y23)).multiply((d2.negate().add(d3.add(t23))).sin())).subtract(q2);

		
		DerivativeStructure[] equations = Bus.createEqauations(buses);
		System.out.println(equations[0].getValue() + ">>>>>>" + fp2.getValue());
		System.out.println(equations[2].getValue() + ">>>>>>" + fp3.getValue());
		System.out.println(equations[1].getValue() + ">>>>>>" + fq2.getValue());
		System.out.println(equations[2].getValue() + ">>>>>>" + fp3.getValue());
		System.out.println(equations[1].getValue() + ">>>>>>" + fq2.getValue());
		
		System.out.println(equations[0].getPartialDerivative(0, 0, 1, 0, 0, 0)+ ">>>>>>" + fp2.getPartialDerivative(1, 0, 0, 0,0,0));
		System.out.println(equations[0].getPartialDerivative(0, 0, 0, 0, 1, 0)+ ">>>>>>" + fp2.getPartialDerivative(0, 1, 0, 0,0,0));
		System.out.println(equations[0].getPartialDerivative(0, 0, 0, 1, 0, 0)+ ">>>>>>" + fp2.getPartialDerivative(0, 0, 1, 0,0,0));
		System.out.println(equations[1].getPartialDerivative(0, 0, 1, 0, 0, 0)+ ">>>>>>" + fq2.getPartialDerivative(1, 0, 0, 0,0,0));
		System.out.println(equations[1].getPartialDerivative(0, 0, 0, 0, 1, 0)+ ">>>>>>" + fq2.getPartialDerivative(0, 1, 0, 0,0,0));
		System.out.println(equations[1].getPartialDerivative(0, 0, 0, 1, 0, 0)+ ">>>>>>" + fq2.getPartialDerivative(0, 0, 1, 0,0,0));

		for (int i = 1; i <= 5; i++) {
			equations = Bus.createEqauations(buses);
			double[] unknowns0 = { buses[1].delta.getValue(), buses[2].delta.getValue(), buses[1].voltage.getValue() };
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
			System.out.println(fx0);
			RealMatrix JI = MatrixUtils.inverse(J);
			X0 = X0.subtract(JI.multiply(fx0));
			buses[1].delta = new DerivativeStructure(params, order, 2, X0.getEntry(0, 0));
			buses[1].voltage = new DerivativeStructure(params, order, 3, X0.getEntry(2, 0));
			buses[2].delta = new DerivativeStructure(params, order, 4, X0.getEntry(1, 0));
			System.out.printf("\td2 = %7.6f \t iteration %d %n", X0.getEntry(0, 0), i);
			System.out.printf("\td3 = %7.6f \t iteration %d %n", X0.getEntry(1, 0), i);
			System.out.printf("\tv2 = %7.6f \t iteration %d %n", X0.getEntry(2, 0), i);
			System.out.printf("%s%n", "--------------------------------------------");
		}

//		for (int i = 1; i <= 10; i++) {
//			DerivativeStructure d2 = new DerivativeStructure(params, order, 0, d2RealValue);
//			DerivativeStructure d3 = new DerivativeStructure(params, order, 1, d3RealValue);
//			DerivativeStructure v2 = new DerivativeStructure(params, order, 2, v2RealValue);
//			DerivativeStructure v1 = new DerivativeStructure(params, order, 3, 1.0);
//			DerivativeStructure d1 = new DerivativeStructure(params, order, 3, 0);
//			DerivativeStructure v3 = new DerivativeStructure(params, order, 3, 1.01);
//
//			DerivativeStructure fp2 = v2.multiply(v1.multiply(Y21)).multiply(d2.negate().add(d1.add(t21)).cos())
//					.add(v2.pow(2).multiply(Y22 * FastMath.cos(t22)))
//					.add(v2.multiply(v3.multiply(Y23)).multiply((d2.negate().add(d3.add(t23))).cos())).subtract(p2);
//			DerivativeStructure fp3 = (d3.negate().add(d1.add(t31)).cos()).multiply(v1.multiply(v3.multiply(Y31)))
//					.add((d3.negate().add(d2).add(t32).cos()).multiply(v2.multiply(v3.multiply(Y32))))
//					.add(v3.pow(2).multiply(Y33).multiply(Math.cos(t33))).subtract(p3);
//			DerivativeStructure fq2 = v2.multiply(v1.multiply(-Y21)).multiply((d2.negate().add(d1.add(t21))).sin())
//					.add(v2.pow(2).multiply(-Y22).multiply(FastMath.sin(t22)))
//					.add(v3.multiply(v2.multiply(-Y23)).multiply((d2.negate().add(d3.add(t23))).sin())).subtract(q2);
//
//			double[] unknowns = { d2.getValue(), d3.getValue(), v2.getValue() };
//			X = new Array2DRowRealMatrix(unknowns);
//			double[][] jacobs = {
//					{ fp2.getPartialDerivative(1, 0, 0, 0), fp2.getPartialDerivative(0, 1, 0, 0),
//							fp2.getPartialDerivative(0, 0, 1, 0) },
//					{ fp3.getPartialDerivative(1, 0, 0, 0), fp3.getPartialDerivative(0, 1, 0, 0),
//							fp3.getPartialDerivative(0, 0, 1, 0) },
//					{ fq2.getPartialDerivative(1, 0, 0, 0), fq2.getPartialDerivative(0, 1, 0, 0),
//							fq2.getPartialDerivative(0, 0, 1, 0) } };
//			RealMatrix J = new Array2DRowRealMatrix(jacobs);
//			double[] functions = { fp2.getValue(), fp3.getValue(), fq2.getValue() };
//			RealMatrix fx = new Array2DRowRealMatrix(functions);
//
//			// System.out.println(J);
//			// LUDecomposition lu =new LUDecomposition(J);
//			RealMatrix JI = MatrixUtils.inverse(J);
//			// System.out.println(lu.getDeterminant());
//
//			X = X.subtract(JI.multiply(fx));
//
//			d2RealValue = X.getEntry(0, 0);
//			d3RealValue = X.getEntry(1, 0);
//			v2RealValue = X.getEntry(2, 0);
//
//			System.out.printf("\td2 = %7.6f \t iteration %d %n", X.getEntry(0, 0), i);
//			System.out.printf("\td3 = %7.6f \t iteration %d %n", X.getEntry(1, 0), i);
//			System.out.printf("\tv2 = %7.6f \t iteration %d %n", X.getEntry(2, 0), i);
//			System.out.printf("%s%n", "--------------------------------------------");
//
//		}

//		for (int i = 1; i <= 10; i++) {
//		DerivativeStructure d2 = new DerivativeStructure(params, order, 0, d2RealValue);
//		DerivativeStructure d3 = new DerivativeStructure(params, order, 1, d3RealValue);
//		DerivativeStructure v2 = new DerivativeStructure(params, order, 2, v2RealValue);
//		DerivativeStructure v1 = new DerivativeStructure(params, order, 3, 1.0);
//		DerivativeStructure d1 = new DerivativeStructure(params, order, 4, 0);
//		DerivativeStructure v3 = new DerivativeStructure(params, order, 5, 1.01);
//
//		DerivativeStructure fp2 = v2.multiply(v1.multiply(Y21)).multiply(d2.negate().add(d1.add(t21)).cos())
//				.add(v2.pow(2).multiply(Y22 * FastMath.cos(t22)))
//				.add(v2.multiply(v3.multiply(Y23)).multiply((d2.negate().add(d3.add(t23))).cos())).subtract(p2);
//		DerivativeStructure fp3 = (d3.negate().add(d1.add(t31)).cos()).multiply(v1.multiply(v3.multiply(Y31)))
//				.add((d3.negate().add(d2).add(t32).cos()).multiply(v2.multiply(v3.multiply(Y32))))
//				.add(v3.pow(2).multiply(Y33).multiply(Math.cos(t33))).subtract(p3);
//		DerivativeStructure fq2 = v2.multiply(v1.multiply(-Y21)).multiply((d2.negate().add(d1.add(t21))).sin())
//				.add(v2.pow(2).multiply(-Y22).multiply(FastMath.sin(t22)))
//				.add(v3.multiply(v2.multiply(-Y23)).multiply((d2.negate().add(d3.add(t23))).sin())).subtract(q2);
//
//		double[] unknowns = { d2.getValue(), d3.getValue(), v2.getValue() };
//		X = new Array2DRowRealMatrix(unknowns);
//		double[][] jacobs = {
//				{ fp2.getPartialDerivative(1, 0, 0, 0, 0, 0), fp2.getPartialDerivative(0, 1, 0, 0, 0, 0),
//						fp2.getPartialDerivative(0, 0, 1, 0, 0, 0) },
//				{ fp3.getPartialDerivative(1, 0, 0, 0, 0, 0), fp3.getPartialDerivative(0, 1, 0, 0, 0, 0),
//						fp3.getPartialDerivative(0, 0, 1, 0, 0, 0) },
//				{ fq2.getPartialDerivative(1, 0, 0, 0, 0, 0), fq2.getPartialDerivative(0, 1, 0, 0, 0, 0),
//						fq2.getPartialDerivative(0, 0, 1, 0, 0, 0) } };
//		RealMatrix J = new Array2DRowRealMatrix(jacobs);
//		double[] functions = { fp2.getValue(), fp3.getValue(), fq2.getValue() };
//		RealMatrix fx = new Array2DRowRealMatrix(functions);
//
//		// System.out.println(J);
//		// LUDecomposition lu =new LUDecomposition(J);
//		RealMatrix JI = MatrixUtils.inverse(J);
//		// System.out.println(lu.getDeterminant());
//
//		X = X.subtract(JI.multiply(fx));
//
//		d2RealValue = X.getEntry(0, 0);
//		d3RealValue = X.getEntry(1, 0);
//		v2RealValue = X.getEntry(2, 0);
//
//		System.out.printf("\td2 = %7.6f \t iteration %d %n", X.getEntry(0, 0), i);
//		System.out.printf("\td3 = %7.6f \t iteration %d %n", X.getEntry(1, 0), i);
//		System.out.printf("\tv2 = %7.6f \t iteration %d %n", X.getEntry(2, 0), i);
//		System.out.printf("%s%n", "--------------------------------------------");
//
//	}

	}
}