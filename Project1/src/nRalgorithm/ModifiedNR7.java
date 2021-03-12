package nRalgorithm;

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;

public class ModifiedNR7 {
	private static Random random = new Random();

	public static void main(String[] args) {

		int nPar = 100;
		Particle particles[] = new Particle[nPar];
		for (int i = 0; i < nPar; i++) {
			particles[i] = new Particle(
					new double[] { 200 * random.nextDouble() - 100, 200 * random.nextDouble() - 100 },
					new double[] { 200 * random.nextDouble() - 100, 200 * random.nextDouble() - 100 });

			particles[i].sol = new Solution2();

		}

		ToyFitnessFunction tff = new ToyFitnessFunction();


		int iteration = 1;
		int max_iter = 1000;
		double bestFitness = Double.MAX_VALUE;
		double[] bestPosition = new double[2];
		int n0 = 10;
	
		
		for (Particle p : particles)
			p.setFitness(tff.getFitness(p, n0));

		for (iteration = 1; iteration < max_iter; iteration++) {
			boolean flag = false;
			for (Particle particle : particles) {
				double[] pOldPos = particle.getPosition().clone();

//				particle.setFitness(tff.getFitness(particle, (int) 100+(iteration%20)*10/particles.length));
				
				particle.setFitness(tff.getFitness(particle, 20));
				ModifiedNR7.makeAdditionalReps(particles, 100+(iteration%20)*10);

				if (particle.getFitness() <= particle.getBestFitness()) {
					particle.setBestFitness(particle.getFitness());
					particle.setBestPosition(pOldPos);

					if (particle.getFitness() <= bestFitness) {
						flag = true;
						bestFitness = particle.getFitness();
						bestPosition = pOldPos;

					}
				}
				double[] position = particle.getPosition();
				double[] speed = particle.getSpeed();
				for (int m = 0; m < position.length; m++) {
					position[m] += speed[m];

				}

				// Updates the particle speed.
				for (int m = 0; m < speed.length; m++) {
					speed[m] = getNewParticleSpeedForIndex(particle, m, iteration, max_iter, bestPosition);
				}

			}
			if (flag) {
				printBest(bestPosition, iteration);
				System.out.printf("%.4f\n",bestFitness);
			} else {
				if (iteration % 100 != 0)
					System.out.print(".");
				else
					System.out.println();
			}
		}
		System.out.println("\nEND");

	}

	private static double getNewParticleSpeedForIndex(Particle particle, int index, int i, int max_iter,
			double[] bestPosition) {

		double c1 = 2.05;
		double c2 = 2.05;

		double x = (((double) max_iter + 1 - i) / ((double) max_iter + 1))
				* (2 / Math.abs((2 - c1 - c2 - Math.sqrt(Math.pow(c1 + c2, 2) - 4 * (c1 + c2)))));

		return x * (particle.getSpeed()[index]
				+ (randomizePercentage(c1) * (particle.getBestPosition()[index] - particle.getPosition()[index]))
				+ (randomizePercentage(c2) * (bestPosition[index] - particle.getPosition()[index])));

	}

	private static double randomizePercentage(double value) {
		return random.nextDouble() * value;
	}

	private static void printBest(double[] best, int f) {
		System.out.printf("\nBest position found at iteration : %d \t", f);
		for (int i = 0; i < best.length; i++) {
			System.out.printf("%.6f\t", best[i]);
		}
	}

	private static int findBestInd(Particle[] particles) {

		double min = Double.MAX_VALUE;
		int bestInd = 0;

		for (int i = 0; i < particles.length; i++) {
			Solution2 temp = particles[i].sol;
			double mean = temp.getMean();
			if (mean < min) {
				min = mean;
				bestInd = i;

			}

		}

		return bestInd;
	}

	public static int[] findRep(Particle[] particles, int budget) {
		double[] ratio = new double[particles.length];
		int bestInd = findBestInd(particles);
		double bestMean = particles[bestInd].sol.getMean();
		int[] adReps = new int[particles.length];

		for (int i = 0; i < particles.length; i++) {
			if (i != bestInd) {
				Solution2 temp = particles[i].sol;
				double mean = temp.getMean();
				double sd = temp.getSD();
				ratio[i] = Math.pow(sd, 2) / Math.pow(mean - bestMean, 2);

			}
		}
		double tempSum = 0;
		for (int i = 0; i < particles.length; i++) {
			if (i != bestInd) {
				tempSum = Math.pow(ratio[i], 2) / Math.pow(particles[i].sol.getSD(), 2);

			}
		}
		ratio[bestInd] = particles[bestInd].sol.getSD() * Math.sqrt(tempSum);

		double sumRatios = 0;

		for (double i : ratio)
			sumRatios += i;

		for (int i = 0; i < ratio.length; i++)
			ratio[i] = ratio[i] / sumRatios;

		for (int i = 0; i < particles.length; i++) {
			adReps[i] = (int) Math.round(budget * ratio[i] > 1 ? budget * ratio[i] : 0);
		}

		return adReps;
	}

	public static void makeAdditionalReps(Particle[] particles, int budget) {
		int[] adReps = ModifiedNR7.findRep(particles, budget);
		NormalDistribution normal = new NormalDistribution(0, 5.0);

		for (int i = 0; i < particles.length; i++) {
			double x1 = particles[i].getPosition()[0];
			double x2 = particles[i].getPosition()[1];
			int rep = adReps[i];
			particles[i].sol.rep += rep;

			if (rep != 0) {
				for (int j = 0; j < rep; j++) {
					double value = (4 * Math.pow(x1, 2) - 2.1 * Math.pow(x1, 4) + Math.pow(x1, 6) / 3 + x1 * x2
							- 4 * Math.pow(x2, 2) + 4 * Math.pow(x2, 4) + 1.0316) + normal.sample();

					particles[i].sol.fitness.add(value);

				}

				particles[i].setFitness(particles[i].sol.getMean());
			}
		}
	}

}
