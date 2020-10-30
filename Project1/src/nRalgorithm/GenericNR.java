package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;

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

}
