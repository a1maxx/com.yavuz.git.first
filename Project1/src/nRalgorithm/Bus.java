package nRalgorithm;

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class Bus {

	public int type;
	public RealMatrix admittance;
	public RealMatrix theta;
	double p, q,nominal_p,nominal_q;
	public DerivativeStructure voltage;
	public DerivativeStructure delta;
	public int index;
	public int order = 2;
	Random random = new Random();
//	public double initialValue_d = random.nextDouble();
//	public double initialValue_v = random.nextDouble()*2;
//	public double initialValue_d = 0.45;
//	public double initialValue_v = 1.04;
	public double initialValue_d = 0.0;
	public double initialValue_v = 1.0;
	boolean slack = false;
	Complex cVolt;
	double mp;
	double nq;
	double qMax=0.0;

	// Constructor of the bus objects
	Bus(int type, double[][] admittance, double[][] theta, int index, double P, double Q) {
		this.type = type;
		
		//This needs to be generic
		int params = 6;
		
		this.admittance = new Array2DRowRealMatrix(admittance);
		this.theta = new Array2DRowRealMatrix(theta);
		this.index = index;
		this.p = P;
		this.q = Q;
		this.mp = 0.5;
		this.nq = 0.5;
		this.delta = new DerivativeStructure(params, order, index, initialValue_d);
		this.voltage = new DerivativeStructure(params, order, index + 1, initialValue_v);
		this.cVolt = Admittance.polarToComplex(voltage.getValue(), delta.getValue());

	}

	Bus(int type, double[][] admittance, double[][] theta, int index, double P, double Q,double mp,double nq) {
		this.type = type;
		int params = 6;
		this.admittance = new Array2DRowRealMatrix(admittance);
		this.theta = new Array2DRowRealMatrix(theta);
		this.index = index;
		this.p = P;
		this.q = Q;
		this.mp = mp;
		this.nq = nq;
		this.delta = new DerivativeStructure(params, order, index, initialValue_d);
		this.voltage = new DerivativeStructure(params, order, index + 1, initialValue_v);
		this.cVolt = Admittance.polarToComplex(voltage.getValue(), delta.getValue());

	}


	// Dummy Constructor
	public Bus() {

	}

	// Creation of the power flow equations for a given array containing Bus objects
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
				equations[l] = equationQ.subtract(buses[k].q);
				;
				l++;
			}

		}
		return equations;
	}

	// Determination of the unknowns
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
	// toString method for class Bus
	public String toString() {

		StringBuilder sb = new StringBuilder();
		sb.append("This is Bus with the type of ");
		sb.append(this.type);
		sb.append("\nIndex of the bus is ");
		sb.append(this.index);
		return sb.toString();

	}
	// Function that generates P and Q functions
	static DerivativeStructure[] createEquations2(ArrayList<Bus> buses) {
		int l = 0;
		int size=0;
		for(Bus b: buses) {
			if(b.type==0) {
				size = size +2;
			}else if(b.type==1)
				size++;
		}
		DerivativeStructure equations[] = new DerivativeStructure[size];
		for (int k = 0; k < buses.size(); k++) {
			DerivativeStructure equationP = new DerivativeStructure(6, 2);
			DerivativeStructure equationQ = new DerivativeStructure(6, 2);
			for (int i = 0; i < buses.size(); i++) {
				if (buses.get(k).admittance.getEntry(k, i) != 0) {
					if (buses.get(k).type == 0) {
						equationP = equationP.add(buses.get(k).voltage.multiply(buses.get(k).admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
										.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).cos()));

						equationQ = equationQ.add(buses.get(k).voltage.multiply(-buses.get(k).admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
										.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).sin()));

					} else if (buses.get(k).type == 1) {
						equationP = equationP.add(buses.get(k).voltage.multiply(buses.get(k).admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
										.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).cos()));

					} else if (buses.get(k).type == 2) {
						equationQ = equationQ.add(buses.get(k).voltage.multiply(buses.get(k).admittance.getEntry(k, i))
								.multiply(buses.get(i).voltage).multiply((buses.get(k).delta.negate()
										.add(buses.get(k).theta.getEntry(k, i)).add(buses.get(i).delta)).sin()));

					} else if (buses.get(k).type == 3) {
						buses.get(k).slack = true;
					}

				}
			}
			if (buses.get(k).type == 0) {
				equations[l] = equationP.subtract(buses.get(k).p);
				l++;
				equations[l] = equationQ.subtract(buses.get(k).q);
				l++;
			} else if (buses.get(k).type == 1) {
				equations[l] = equationP.subtract(buses.get(k).p);
				l++;
			} else if (buses.get(k).type == 2) {
				equations[l] = equationQ.subtract(buses.get(k).q);
				;
				l++;
			}

		}
		return equations;
	}

	Bus(int type,int index,int params, double P, double Q,double mp,double nq,double qmax) {
		this.type = type;
		this.index = index;
		this.p = P;
		this.q = Q;
		this.mp = mp;
		this.nq = nq;
		this.delta = new DerivativeStructure(params, order, index, initialValue_d);
		this.voltage = new DerivativeStructure(params, order, index + 1, initialValue_v);
		this.cVolt = Admittance.polarToComplex(voltage.getValue(), delta.getValue());
		this.qMax = qmax;
	}
	
	Bus(int type,int index,int params, double P, double Q,double V) {
		this.type = type;
		this.index = index;
		this.p = P;
		this.q = Q;
		this.delta = new DerivativeStructure(params, order, index, initialValue_d);
		this.voltage = new DerivativeStructure(params, order, index + 1, V);
		this.cVolt = Admittance.polarToComplex(voltage.getValue(), delta.getValue());

	}
	
	Bus(int type,int index,int params, double P, double Q) {
		this.type = type;
		this.index = index;
		this.nominal_p=P;
		this.nominal_q=Q;
		this.p = P;
		this.q = Q;
		this.delta = new DerivativeStructure(params, order, index, initialValue_d);
		this.voltage = new DerivativeStructure(params, order, index + 1, initialValue_v);
		this.cVolt = Admittance.polarToComplex(voltage.getValue(), delta.getValue());

	}
	Bus(int type,int index,int params, double P) {
		double PF = 0.328;
		this.type = type;
		this.index = index;
		this.nominal_p=P;
		this.nominal_q=P*PF;;
		this.p = P;
		this.q = P*PF;
		this.delta = new DerivativeStructure(params, order, index, initialValue_d);
		this.voltage = new DerivativeStructure(params, order, index + 1, initialValue_v);
		this.cVolt = Admittance.polarToComplex(voltage.getValue(), delta.getValue());

	}

	
}