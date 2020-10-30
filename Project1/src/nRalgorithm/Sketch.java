package nRalgorithm;

import java.util.ArrayList;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class Sketch {

	public static void main(String[] args) {
		ArrayList<Double> unknowns = new ArrayList<Double>();
		if (true) {
			Double temp = 1.0;
			unknowns.add(temp);
			temp = 2.0;
			unknowns.add(temp);
			System.out.println(unknowns);
		}
		RealMatrix X0 = new Array2DRowRealMatrix();
		X0 = new Array2DRowRealMatrix(GenericNR.convertArray(unknowns));
		System.out.println(X0);
	}

}
