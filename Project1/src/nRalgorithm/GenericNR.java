package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.linear.RealMatrix;

public class GenericNR {

	public static void main(String[] Args) {

	}

	static ArrayList<Integer[]> createOrders(ArrayList<Bus> buses) {

		ArrayList<Integer[]> orders = new ArrayList<Integer[]>();

		for (Bus b : buses) {

			switch (b.type) {
			case 0:
				Integer[] temp = new Integer[buses.size() * 2];
				Arrays.fill(temp, 0);
				temp[b.index] = 1;
				orders.add(temp);
				temp = new Integer[buses.size() * 2];
				Arrays.fill(temp, 0);
				temp[b.index + 1] = 1;
				orders.add(temp);
				break;
			case 1:
				temp = new Integer[buses.size() * 2];
				Arrays.fill(temp, 0);
				temp[b.index] = 1;
				orders.add(temp);
				break;
			case 2:
				temp = new Integer[buses.size() * 2];
				Arrays.fill(temp, 0);
				temp[b.index + 1] = 1;
				orders.add(temp);
				break;
			default:

			}

		}
		return orders;
	}

	static ArrayList<Double> calculateUnknows(ArrayList<Bus> buses) {
		ArrayList<Double> unknowns = new ArrayList<Double>();
		for (Bus b : buses) {

			switch (b.type) {
			case 0:
				Double temp = b.delta.getValue();
				unknowns.add(temp);
				temp = b.voltage.getValue();
				unknowns.add(temp);
				break;
			case 1:
				temp = b.delta.getValue();
				unknowns.add(temp);
				break;
			case 2:
				temp = b.voltage.getValue();
				unknowns.add(temp);
				break;
			default:

			}

		}

		return unknowns;
	}

	static double[] convertArray(ArrayList<Double> unknowns) {
		double[] dArray = new double[unknowns.size()];

		for (int i = 0; i < unknowns.size(); i++) {
			dArray[i] = unknowns.get(i);
		}

		return dArray;
	}

	static double[][] calculateJacobian(ArrayList<Integer[]> orders, DerivativeStructure[] equations) {
		double[][] jacobian = new double[equations.length][orders.size()];
		for (int i = 0; i < equations.length; i++) {
			for (int j = 0; j < orders.size(); j++) {
				jacobian[i][j] = equations[i].getPartialDerivative(ArrayUtils.toPrimitive(orders.get(j)));
			}
		}
		return jacobian;
	}

	static double[] calculateFunctions(DerivativeStructure[] equations) {
		double[] values = new double[equations.length];
		for (int i = 0; i < equations.length; i++) {
			values[i] = equations[i].getValue();
		}
		return values;
	}

	static void updateUnknowns(RealMatrix X0, int params, int order, ArrayList<Integer[]> orders,
			ArrayList<Bus> buses) {
		int j = 0;
		for (Integer[] o : orders) {
			int k = Arrays.asList(o).indexOf(1);
			if (k % 2 == 0) {
				buses.get(k / 2).delta = new DerivativeStructure(params, order, k, X0.getEntry(j, 0));
				j++;
			} else if (k % 2 == 1) {
				buses.get(k / 2).voltage = new DerivativeStructure(params, order, k, X0.getEntry(j, 0));
				j++;
			}

		}
	}
	

}
