package nRalgorithm;

public class Solution {

	public int rep;
	public double[] fitness;
	public double mean;
	public double sd;

	public Solution(double[] fitness, int rep) {
		this.rep = rep;
		this.mean = getMean();
		this.sd = getSD();
		this.fitness = fitness;

	}

	public  double getSD() {

		double mean = this.getMean();
		double sum = 0;
		for (int i = 0; i < this.fitness.length; i++) {
			sum += (this.fitness[i] - mean) / (this.fitness.length - 1);
		}

		return Math.sqrt(sum);

	}

	public  double getMean() {		
		double sum = 0;
		for (int i = 0; i < this.fitness.length; i++) {
			sum += this.fitness[i];
		}

		return sum / this.fitness.length;
	}

	public static int findMaxIndex(double[] arr) {

		int maxInd = 0;
		double max = 0;
		for (int i = 0; i < arr.length; i++) {
			if (arr[i] > max) {
				max = arr[i];
				maxInd = i;
			}
		}

		return maxInd;
	}

	public static int findMinIndex(double[] arr) {
		int minInd = 0;
		double min = Double.MAX_VALUE;
		for (int i = 0; i < arr.length; i++) {
			if (arr[i] < min) {
				min = arr[i];
				minInd = i;
			}
		}

		return minInd;
	}

}
