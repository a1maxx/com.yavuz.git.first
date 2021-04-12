package nRalgorithm;

import java.util.ArrayList;

import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.FastMath;

public class ModifiedNR7 {
	private static Random random = new Random();

	public static void main(String[] args) {

		int nPar = 5;
		Particle particles[] = new Particle[nPar];
		for (int i = 0; i < nPar; i++) {
			particles[i] = new Particle(
					new double[] { 200 * random.nextDouble() - 100, 200 * random.nextDouble() - 100 },
					new double[] { 200 * random.nextDouble() - 100, 200 * random.nextDouble() - 100 });

			particles[i].sol = new Solution2();

		}

		ToyFitnessFunction tff = new ToyFitnessFunction();

		int iteration = 1;
		int max_iter = 100;
		double bestFitness = Double.MAX_VALUE;
		double[] bestPosition = new double[2];
		int n0 = 100;

		// Initialization
		for (Particle particle2 : particles) {
			particle2.setFitness(tff.getFitness(particle2, n0));
			particle2.setBestPosition(particle2.getPosition());
			particle2.setBestFitness(particle2.getFitness());
		}

		bestFitness = particles[ModifiedNR7.findBestInd(particles)].getFitness();
		bestPosition = particles[ModifiedNR7.findBestInd(particles)].getPosition().clone();

		// End of initialization

		Particle best = new Particle();
		best.setPosition(bestPosition);
		best.setFitness(bestFitness);
		best.sol = particles[ModifiedNR7.findBestInd(particles)].sol;

		
		
		for (int i=0;i<particles.length;i++) {
			particles[i].setFitness(tff.getFitness(particles[i], n0));
		}
		
//Case1
		// PSO start
		for (iteration = 1; iteration < max_iter; iteration++) {
			boolean flag = false;
			for (Particle particle : particles) {

				
//				particle.setFitness(tff.getFitness(particle, (int) 700/particles.length));

//Case 2
//				Res res = null;
//				ModifiedNR7.makeAdditionalReps(particles, 100, bestFitness, res);

//Case 3
				Res res = ModifiedNR7.findCase(particles, bestFitness);
				ModifiedNR7.makeAdditionalReps2(particles, 100, bestFitness, res, best);

				double[] pOldPos = particle.getPosition().clone();

				if (particle.getFitness() <= particle.getBestFitness()) {
					particle.setBestFitness(particle.getFitness());
					particle.setBestPosition(pOldPos);

					if (particle.getFitness() <= best.getFitness()) {
						flag = true;
						bestFitness = particle.getFitness();
						bestPosition = pOldPos;

						best.setPosition(bestPosition);
						best.setFitness(bestFitness);
						
						best.sol = particles[ModifiedNR7.findBestInd(particles)].sol;

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
				System.out.printf("%.4f\n", bestFitness);
			} else {
				if (iteration % 100 != 0)
					System.out.print(".");
				else
					System.out.println();
			}
		}
		System.out.println("\nEND");

	}

	public static double randomizePercentage(double value) {
		return random.nextDouble() * value;
	}

	public static void printBest(double[] best, int f) {
		System.out.printf("\nBest position found at iteration : %d \t", f);
		for (int i = 0; i < best.length; i++) {
			System.out.printf("%.6f,\t", best[i]);
		}
	}

	public static int findBestInd(Particle[] particles) {

		double min = Double.MAX_VALUE;
		int bestInd = 0;

		for (int i = 0; i < particles.length; i++) {

			double fitness = particles[i].getFitness();
			if (fitness < min) {
				min = fitness;
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
				tempSum += Math.pow(ratio[i], 2) / Math.pow(particles[i].sol.getSD(), 2);

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

	public static void makeAdditionalReps(Particle[] particles, int budget, double currentBest, Res res) {
		int[] adReps = ModifiedNR7.findRep(particles, budget);

//		int[] adReps = ModifiedNR7.findRep2(particles, budget, currentBest, res);

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

	public static void makeAdditionalReps2(Particle[] particles, int budget, double currentBest, Res res,
			Particle best) {
//		int[] adReps = ModifiedNR7.findRep(particles, budget);
		NormalDistribution normal = new NormalDistribution(0, 5.0);
		int[] adReps = ModifiedNR7.findRep2(particles, budget, currentBest, res,best);

		for (int i = 0; i > adReps[particles.length]; i++) {
			double x1 = best.getPosition()[0];
			double x2 = best.getPosition()[1];
			double value = (4 * Math.pow(x1, 2) - 2.1 * Math.pow(x1, 4) + Math.pow(x1, 6) / 3 + x1 * x2
					- 4 * Math.pow(x2, 2) + 4 * Math.pow(x2, 4) + 1.0316) + normal.sample();
			best.sol.fitness.add(value);
		}
		best.setFitness(best.sol.getMean());
		
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
				particles[i].sol.sd = particles[i].sol.getSD();
			}
		}
	}

	public static int[] findRep2(Particle[] particles, int budget, double currentBest, Res res, Particle best) {
		double[] ratio = new double[particles.length + 1];
	
	
		double bestMean = best.sol.getMean();
		int[] adReps = new int[particles.length + 1];

		// Case of not having 'allegedly' better particle than the currentBest
		if (!res.a) {

			for (int i = 0; i < particles.length; i++) {

				Solution2 temp = particles[i].sol;
				double mean = temp.getMean();
				double sd = temp.getSD();

				if (mean <= particles[i].getBestFitness() && I(particles[i]) > G(particles[i], currentBest))
					ratio[i] = Math.pow(sd, 2) / Math.pow(mean - particles[i].getBestFitness(), 2);
				else if (mean <= particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
					ratio[i] = Math.pow(sd, 2) / Math.pow(mean - bestMean, 2);
				else if (mean > particles[i].getBestFitness() && I(particles[i]) > G(particles[i], currentBest))
					ratio[i] = Math.pow(sd, 2) / Math.pow(mean - particles[i].getBestFitness(), 2);
				else if (mean > particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
					ratio[i] =  Math.pow(sd, 2) / Math.pow(mean - bestMean, 2);
				else 
					System.out.println("ERROR");

			}
			
			for(int i=0;i<ratio.length;i++) {
				if(Double.isInfinite(ratio[i]))
					ratio[i]=0.02;
					
			}
			
			double tempSum1 = 0;
			double tempSum2 = 0;
			for (int i = 0; i < particles.length; i++) {
				Solution2 temp = particles[i].sol;
				double mean = temp.getMean();
				double sd = temp.getSD();

				if (mean <= particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
					tempSum1 += Math.pow(ratio[i], 2) / Math.pow(sd, 2);
				else if (mean > particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
					tempSum2 += Math.pow(ratio[i], 2) / Math.pow(sd, 2);

			}
			ratio[particles.length] = best.sol.getSD() * Math.sqrt(tempSum1 + tempSum2);

			// Case of having better particle than the currentBest
		} else if (res.a) {

			int allegedBest = res.i;
			for (int i = 0; i < particles.length; i++) {
				Solution2 temp = particles[i].sol;
				double mean = temp.getMean();
				double sd = temp.getSD();
				if (i != allegedBest) {
					if (mean <= particles[i].getBestFitness() && I(particles[i]) > G(particles[i], currentBest))
						ratio[i] = Math.pow(sd, 2) / Math.pow((mean - particles[i].getBestFitness()), 2);
					else if (mean <= particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
						ratio[i] = Math.pow(sd, 2) / Math.pow(mean - res.fitness, 2);
					else if (mean > particles[i].getBestFitness() && I(particles[i]) > G(particles[i], currentBest))
						ratio[i] = Math.pow(sd, 2) / Math.pow(mean - particles[i].getBestFitness(), 2);
					else if (mean > particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
						ratio[i] = ratio[i] = Math.pow(sd, 2) / Math.pow(mean - res.fitness, 2);
				}
			}
			ratio[particles.length] = Math.pow(best.sol.getSD(), 2)
					/ Math.pow((best.getFitness() - res.fitness), 2);

			double tempSum1 = 0;
			double tempSum2 = 0;
			for (int i = 0; i < particles.length; i++) {
				Solution2 temp = particles[i].sol;
				double mean = temp.getMean();
				double sd = temp.getSD();
				if (i != allegedBest) {
					if (mean <= particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
						tempSum1 += Math.pow(ratio[i], 2) / Math.pow(sd, 2);
					else if (mean > particles[i].getBestFitness() && I(particles[i]) <= G(particles[i], currentBest))
						tempSum2 += Math.pow(ratio[i], 2) / Math.pow(sd, 2);

				}

			}

			double term3 = Math.pow(ratio[particles.length], 2) / Math.pow(best.sol.getSD(), 2);
			ratio[allegedBest] = particles[allegedBest].sol.getSD() * Math.sqrt(tempSum1 + tempSum2 + term3);

		}

		double sumRatios = 0;

		for (double i : ratio)
			sumRatios += i;

		for (int i = 0; i < ratio.length; i++)
			ratio[i] = ratio[i] / sumRatios;

		for (int i = 0; i < particles.length + 1; i++) {
			adReps[i] = (int) Math.round(budget * ratio[i] > 1 ? budget * ratio[i] : 0);
//			System.out.printf("%d \t",adReps[i]);
		}

		return adReps;
	}

	public static double[] getFitnesses(Particle[] particles) {
		double[] fitnesses = new double[particles.length];
		for (int i = 0; i < particles.length; i++) {
			fitnesses[i] = particles[i].getFitness();

		}
		return fitnesses;

	}

	

	public static double G(Particle p, double best) {

		double term1 = FastMath.pow(best - p.getFitness(), 2);
		double term2 = FastMath.pow(p.sol.getSD(), 2);

		return term1 / term2;
	}

	public static double I(Particle p) {
		
		 double term1 = FastMath.pow(p.getBestFitness() - p.getFitness(), 2);
		 double term2 = 2 * FastMath.pow(p.sol.getSD(), 2);

		return term1 / term2;
	}

	public static Res findCase(Particle[] particles, double best) {

		Res res = new Res();
		res.a = false;
		res.i = -1;
		for (int i = 0; i < particles.length; i++) {
			double fitness = particles[i].getFitness();
			if (fitness < best) {
				res.a = true;
				res.i = i;
				res.fitness = fitness;
				return res;
			}

		}

		return res;

	}

	public static double getNewParticleSpeedForIndex(Particle particle, int index, int i, int max_iter,
			double[] bestPosition) {

		double c1 = 2.05;
		double c2 = 2.05;

		double x = (((double) max_iter + 1 - i) / ((double) max_iter + 1))
				* (2 / Math.abs((2 - c1 - c2 - Math.sqrt(Math.pow(c1 + c2, 2) - 4 * (c1 + c2)))));

		return x * (particle.getSpeed()[index]
				+ (randomizePercentage(c1) * (particle.getBestPosition()[index] - particle.getPosition()[index]))
				+ (randomizePercentage(c2) * (bestPosition[index] - particle.getPosition()[index])));

	}

}
