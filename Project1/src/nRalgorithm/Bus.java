package nRalgorithm;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class Bus {
	private int a;
	private int type;
	private RealMatrix admittance;
	private RealMatrix theta;
	double p, q;
	public DerivativeStructure voltage;
	public DerivativeStructure delta;
	int order = 2;
	double initialValue_d = 0.0;
	double initialValue_v = 1.0;
	boolean slack = false;

	Bus(int type, double[][] admittance, double[][] theta, int index, double P, double Q) {
		this.type = type;
		// int params = admittance.length;
		int params = 6;
		this.admittance = new Array2DRowRealMatrix(admittance);
		this.theta = new Array2DRowRealMatrix(theta);
		this.p = P;
		this.q = Q;
		delta = new DerivativeStructure(params, order, index, initialValue_d);
		voltage = new DerivativeStructure(params, order, index + 1, initialValue_v);

	}

	public Bus() {

	}

	static DerivativeStructure[] createEqauations(Bus[] buses) {
		int l = 0;
		DerivativeStructure equations[] = new DerivativeStructure[buses.length];
		for (int k = 0; k < buses.length; k++) {
			DerivativeStructure equationP = new DerivativeStructure(6, 2);
			DerivativeStructure equationQ = new DerivativeStructure(6, 2);
			for (int i = 0; i < buses.length; i++) {
				if (buses[k].admittance.getEntry(k, i) != 0) {
					if (buses[k].type == 0) {
						equationP = equationP.add(buses[k].voltage.multiply(buses[k].admittance.getEntry(k, i))
								.multiply(buses[i].voltage).multiply(
										(buses[k].delta.negate().add(buses[k].theta.getEntry(k, i)).add(buses[i].delta))
												.cos()));

						equationQ = equationQ.add(buses[k].voltage.multiply(-buses[k].admittance.getEntry(k, i))
								.multiply(buses[i].voltage).multiply(
										(buses[k].delta.negate().add(buses[k].theta.getEntry(k, i)).add(buses[i].delta))
												.sin()));

					} else if (buses[k].type == 1) {
						equationP = equationP.add(buses[k].voltage.multiply(buses[k].admittance.getEntry(k, i))
								.multiply(buses[i].voltage).multiply(
										(buses[k].delta.negate().add(buses[k].theta.getEntry(k, i)).add(buses[i].delta))
												.cos()));

					} else if (buses[k].type == 2) {
						equationQ = equationQ.add(buses[k].voltage.multiply(buses[k].admittance.getEntry(k, i))
								.multiply(buses[i].voltage).multiply(
										(buses[k].delta.negate().add(buses[k].theta.getEntry(k, i)).add(buses[i].delta))
												.sin()));

					} else if (buses[k].type == 3) {
						buses[k].slack = true;
					}

				}
			}
			if (buses[k].type == 0) {
				equations[l] = equationP.subtract(buses[k].p);
				l++;
				equations[l] = equationQ.subtract(buses[k].q);
				l++;
			} else if (buses[k].type == 1) {
				equations[l] = equationP.subtract(buses[k].p);
				l++;
			} else if (buses[k].type == 2) {
				equations[l] = equationQ.subtract(buses[k].q);;
				l++;
			}

		}
		return equations;
	}

	DerivativeStructure[] detectUnknowns(Bus[] buses) {
		DerivativeStructure unknowns[] = new DerivativeStructure[this.admittance.getRowDimension()];
		int params = unknowns.length;
		int order = 2;
		double initialValue_d = 0.0;
		double initialValue_v = 1.0;
		for (int i = 1; i <= buses.length; i++) {
			if (buses[i].type == 1) {
				unknowns[i] = new DerivativeStructure(params, order, i, initialValue_d);

			} else if (buses[i].type == 2) {
				unknowns[i] = new DerivativeStructure(params, order, i, initialValue_v);

			}
		}
		return unknowns;
	}

}