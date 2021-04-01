package nRalgorithm;

import java.util.Arrays;

public class Particle {
	
	
	public Solution2 sol;

	
	
	/**
	 * The current position of this particle.
	 */

	private double[] position;
	
	private int rep;

	/**
	 * The speed of this particle.
	 */
	private double[] speed;

	/**
	 * The fitness of this particle for the current position.
	 */
	private double fitness;
	

	/**
	 * The best position found by this particle.
	 */
	private double[] bestPosition;

	/**
	 * The best fitness found by this particle.
	 */
	private double bestFitness = Double.POSITIVE_INFINITY;


	public Particle(double[] initialPosition, double[] initialSpeed) {
		this.position = initialPosition;
		this.speed = initialSpeed;
		rep = 10;
		this.sol= new Solution2();
	}

	public Particle() {
		this.sol= new Solution2();
	}
	public double[] getPosition() {
		return position;
	}


	public double[] getSpeed() {
		return speed;
	}


	public double getFitness() {
		return fitness;
	}


	public double[] getBestPosition() {
		return bestPosition;
	}


	public double getBestFitness() {
		return bestFitness;
	}


	public void setPosition(double[] position) {
		this.position = position;
	}


	public void setSpeed(double[] speed) {
		this.speed = speed;
	}


	public void setFitness(double fitness) {
		this.fitness = fitness;
	}
	
	public void setSol(double [] fitness){
		
		
	}

	public void setBestPosition(double[] bestPosition) {
		this.bestPosition = bestPosition;
	}


	public void setBestFitness(double bestFitness) {
		this.bestFitness = bestFitness;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(bestFitness);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + Arrays.hashCode(bestPosition);
		temp = Double.doubleToLongBits(fitness);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + Arrays.hashCode(position);
		result = prime * result + Arrays.hashCode(speed);
		return result;
	}


	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Particle other = (Particle) obj;
		if (Double.doubleToLongBits(bestFitness) != Double.doubleToLongBits(other.bestFitness))
			return false;
		if (!Arrays.equals(bestPosition, other.bestPosition))
			return false;
		if (Double.doubleToLongBits(fitness) != Double.doubleToLongBits(other.fitness))
			return false;
		if (!Arrays.equals(position, other.position))
			return false;
		if (!Arrays.equals(speed, other.speed))
			return false;
		return true;
	}


	@Override
	public String toString() {
		return "Particle [position=" + Arrays.toString(position) + ", speed=" + Arrays.toString(speed) + ", \nfitness="
				+ fitness + ", bestPosition=" + Arrays.toString(bestPosition) + ", bestFitness=" + bestFitness + "]";
	}


	public int getRep() {
		return rep;
	}


	public void setRep(int rep) {
		this.rep = rep;
	}
}