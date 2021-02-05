package nRalgorithm;

import java.util.ArrayList;

public class Main2 {

	public static void main(String[] args) {
		ModifiedNR4 nr4 = new ModifiedNR4();
		
		double [] mp = new double[] {0.12225, 0.72107, 0.10940};
		double [] nq = new double[] {0.2556, 0.86734, 0.86643} ;

		
		
		ArrayList<Bus> buses = nr4.createBuses(mp,nq);
		nr4.runMNR(buses);

		
		
		
	}

}
