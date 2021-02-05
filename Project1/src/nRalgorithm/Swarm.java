package nRalgorithm;

import java.util.Arrays;
import java.util.Random;

public class Swarm {
	/**
	 * The particles of this swarm.
	 */
	private Particle[] particles;

	/**
	 * The best position found within the particles of this swarm.
	 */
	private double[] bestPosition;

	/**
	 * The best fitness score found within the particles of this swarm.
	 */
	private double bestFitness = Double.POSITIVE_INFINITY;

	/**
	 * A random generator.
	 */
	private Random random = new Random();

	/**
	 * Instantiates a new Swarm.
	 *
	 * @param numParticles
	 *            the number of particles of the swarm
	 */
	public Swarm(int numParticles) {
		particles = new Particle[numParticles];
		for (int i = 0; i < numParticles; i++) {
			double[] initialParticlePosition = { random.nextDouble() , random.nextDouble() ,
					random.nextDouble() , random.nextDouble() , random.nextDouble() ,
					random.nextDouble() };
			double[] initialParticleSpeed = { random.nextDouble() , random.nextDouble() ,
					random.nextDouble(), random.nextDouble() , random.nextDouble() ,
					random.nextDouble() };
			particles[i] = new Particle(initialParticlePosition, initialParticleSpeed);
		}
	}

	public Particle[] getParticles() {
		return particles;
	}

	public double[] getBestPosition() {
		return bestPosition;
	}

	public double getBestFitness() {
		return bestFitness;
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
		result = prime * result + Arrays.hashCode(particles);
		result = prime * result + ((random == null) ? 0 : random.hashCode());
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
		Swarm other = (Swarm) obj;
		if (Double.doubleToLongBits(bestFitness) != Double.doubleToLongBits(other.bestFitness))
			return false;
		if (!Arrays.equals(bestPosition, other.bestPosition))
			return false;
		if (!Arrays.equals(particles, other.particles))
			return false;
		if (random == null) {
			if (other.random != null)
				return false;
		} else if (!random.equals(other.random))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "Swarm [particles=" + Arrays.toString(particles) + ", bestPosition=" + Arrays.toString(bestPosition)
				+ ", bestFitness=" + bestFitness + ", random=" + random + "]";
	}

}
