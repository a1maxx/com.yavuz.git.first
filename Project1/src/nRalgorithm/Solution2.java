package nRalgorithm;

import java.util.ArrayList;

public class Solution2 {

	public int rep;
	public ArrayList<Double> fitness;

	
	public Solution2() {
		this.rep=0;
		this.fitness = new ArrayList<Double>();
	}
	
	public double getMean() {
		double sum = 0;
		for (Double d : this.fitness) {
			sum += d;
		}

		return sum / this.fitness.size();

	}

	public double getSD() {

		double mean = this.getMean();
		double sum = 0;
		for (int i = 0; i < this.fitness.size(); i++) {
			sum += Math.pow(this.fitness.get(i) - mean, 2);
		}

		return Math.sqrt(sum / this.fitness.size() - 1);

	}

	public static int findMaxIndex(ArrayList<Double> arr) {

		int maxInd = 0;
		double max = 0;
		for (int i = 0; i < arr.size(); i++) {
			if (arr.get(i) > max) {
				max = arr.get(i);
				maxInd = i;
			}
		}

		return maxInd;
	}

	public static int findMinIndex(ArrayList<Double> arr) {
		int minInd = 0;
		double min = Double.MAX_VALUE;
		for (int i = 0; i < arr.size(); i++) {
			if (arr.get(i) < min) {
				min = arr.get(i);
				minInd = i;
			}
		}

		return minInd;
	}

}
