package nRalgorithm;

import java.util.ArrayList;
import java.util.Collections;

public class Solution2 {

	public int rep;
	public ArrayList<Double> fitness;
	public double mean;
	public double sd;
	
	public Solution2() {
		this.rep=0;
		this.fitness = new ArrayList<Double>();
		
	}
	
	public double getMean() {
		double sum = 0;
		for (Double d : this.fitness) {
			sum += d;
		}
		
		this.mean = sum / this.fitness.size();
		return this.mean;

	}

	public double getSD() {

		double mean = this.getMean();
		double sum = 0;
		for (int i = 0; i < this.fitness.size(); i++) {
			sum += Math.pow(this.fitness.get(i) - mean, 2);
		}

		this.sd = Math.sqrt(sum / (this.fitness.size() - 1));

		return this.sd;
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
