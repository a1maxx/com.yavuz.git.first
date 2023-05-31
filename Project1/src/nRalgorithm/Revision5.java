package nRalgorithm;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;

public class Revision5 {

	public static final double loadMean = 0.15;
	public static final double loadSD = loadMean * 0.4;

	double[][] xadmittances;
	double[][] radmittances;
	double[][] theta;
	WeibullDistribution wb = new WeibullDistribution(7.5, 3.5);

	public Revision5() {
		int nofB = 6;
		double[][] xadmittances = new double[nofB][nofB];
		double[][] radmittances = new double[nofB][nofB];
		File file = new File("case6ww.txt");
		try {
			Scanner scan = new Scanner(file);
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					xadmittances[i][j] = scan.nextDouble();

				}
			}
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					radmittances[i][j] = scan.nextDouble();
				}
			}
			scan.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		this.xadmittances = xadmittances;
		this.radmittances = radmittances;

	}

	@SuppressWarnings("unchecked")
	public static void main(String[] args) {

		Revision5 mnr8 = new Revision5();

		int npar = 0; // Number of particles
		int max_iter = 0; // Maximum Iterations
		int n0 = 0; // Initial replication
		int budget = 0; // budget to allocate
//
//		int[] npars = { 5, 10, 15, 20, 25, 30, 35, 40 };
//		int[] max_iters = { 10, 100, 150, 200, 250, 300, 350, 700 };
//		int[] n0s = { 5, 10, 15, 20, 25, 30, 35, 40 };
//		int[] budgets = { 10, 50, 60, 70, 80, 90, 100, 200 };

//		int[] npars = { 5, 10, 15, 30, 40 };
//		int[] max_iters = { 10, 100, 175, 200};
//		int[] n0s = { 5, 10, 25, 35, 40};
//		int[] budgets = { 10, 25, 60, 90};
		
//		int[] npars = {30, 40, 50 };
//		int[] max_iters = { 50, 100, 125};
//		int[] n0s = {10, 15, 20};
//		int[] budgets = {10, 20, 35};
		
		int[] npars = {50};
		int[] max_iters = {100};
		int[] n0s = {10};
		int[] budgets = {10,100,200,400,800,1600};
		
		double[] bestPosition = new double[4];

		for (int a = 0; a < npars.length; a++) {
			for (int b = 0; b < max_iters.length; b++) {
			
				for (int c = 0; c < n0s.length; c++) {
					
					for (int d = 0; d < budgets.length; d++) {
						n0 = n0s[c]; // Initial replication
						max_iter = max_iters[b]; // Maximum Iterations
						budget = budgets[d]; // budget to allocate
						npar = npars[a]; // Number of particles
						
						Particle[] particles = mnr8.initializeParticles(npar);

						int iteration = 1;
						double bestFitness = Double.MAX_VALUE;

						bestFitness = particles[ModifiedNR7.findBestInd(particles)].getFitness();
						bestPosition = particles[ModifiedNR7.findBestInd(particles)].getPosition().clone();

						particles = mnr8.initializeParticles(npar);
						ArrayList<Double> bestFitnesses = (ArrayList<Double>) particles[ModifiedNR7
								.findBestInd(particles)].sol.fitness.clone();

						Particle best0 = new Particle();
						long start = System.currentTimeMillis();

						for (iteration = 1; iteration < max_iter; iteration++) {
							boolean flag = false;
							
							Particle best = new Particle();
							best.setPosition(bestPosition);
							best.setFitness(bestFitness);
							best.sol.fitness = bestFitnesses;
							mnr8.addSamples(best, n0);
							bestFitness = best.getFitness();

							for (Particle particle2 : particles)
								mnr8.experimentParticle(particle2, n0);

							Res res = ModifiedNR7.findCase(particles, bestFitness);
							mnr8.makeAdditionalReps3(particles, budget, bestFitness, res, best);
							bestFitness = best.getFitness();

							for (Particle particle : particles) {

								double[] pOldPos = particle.getPosition().clone();

								if (particle.getFitness() <= particle.getBestFitness()) {
									particle.setBestFitness(particle.getFitness());
									particle.setBestPosition(pOldPos);

									if (particle.getFitness() <= bestFitness) {
										flag = true;
										bestFitness = particle.getFitness();
										bestFitnesses = (ArrayList<Double>) particle.sol.fitness.clone();
										bestPosition = pOldPos;
									}
								}

								double[] position = particle.getPosition();
								double[] speed = particle.getSpeed();
								
								for (int m = 0; m < position.length; m++) {
									position[m] += speed[m];

								}

								// Update the particle speed

								for (int m = 0; m < speed.length; m++) {
									speed[m] = ModifiedNR7.getNewParticleSpeedForIndex(particle, m, iteration, max_iter,
											bestPosition);
								}

							}

							// If the best particle has changed then print
							if (flag || (iteration+1)==max_iter) {
								best0.setPosition(bestPosition);
								best0.setFitness(bestFitness);
								best0.sol.fitness = bestFitnesses;

//								ModifiedNR7.printBest(bestPosition, iteration);
//								System.out.printf("%.4f\t%.4f\t%d", best.sol.getMean(), best.sol.getSD(), best.sol.fitness.size());
							}
						}
						long end = System.currentTimeMillis();
						ModifiedNR7.printBest(bestPosition, iteration);
//						System.out.printf("%.4f\t%.4f\t%d\n", best0.sol.getMean(), best0.sol.getSD(),
//								best0.sol.fitness.size());
						System.out.printf(
								"\n\n Time elapsed(s) = %d \t npar = %d\t max_iter = %d\t n0 = %d\t budget = %d\t fitness = %.5f \t sd = %.5f\t size= %d\n",
								(end - start)/1000, npar, max_iter, n0, budget, best0.sol.getMean(), best0.sol.getSD(),
								best0.sol.fitness.size());

					}
				}
			}
		}

	}

	public void makeAdditionalReps3(Particle[] particles, int budget, double currentBest, Res res, Particle best) {

		int[] adReps = ModifiedNR7.findRep2(particles, budget, currentBest, res, best);

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;

		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadSD);

		if (adReps[particles.length] > 0)
			this.addSamples(best, adReps[particles.length]);

		for (int f = 0; f < particles.length; f++) {
			int rep = adReps[f];

			if (rep > 0) {
				for (int k = 0; k < rep; k++) {
					int params = 12;
					int order = 2;

					X1 = null;
					fx0 = null;
					wi = 1.0;

					prev_Mismatches = Double.MAX_VALUE;
					flag = true;

					Jacobian = null;

					double[] position = particles[f].getPosition();
					boolean flag2 = true;
					while (flag2) {
						buses = new ArrayList<Bus>();
						buses.add(new Bus(1, 0, params, Math.min(-normal.sample(), 0)));
						buses.add(new Bus(1, 2, params, Math.min(-normal.sample(), 0)));
						buses.add(new Bus(1, 4, params, Math.min(-normal.sample(), 0)));
						buses.add(new Bus(2, 6, params, 0, 0, position[0], position[2], 0.6, 0.8));
						buses.add(new Bus(1, 8, params, generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
						buses.add(new Bus(2, 10, params, 0, 0, position[1], position[3], 0.9, 1.2));
						flag2 = ModifiedNR.checkAdequacy(buses, 2.5);
					}

					ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
					ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

					innerloop: for (int i = 1; i <= 200 && flag; i++) {
						double w0 = 1.00;
						double v0 = 1.01;

						Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances,
								xadmittances, wi);

						ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

						PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

						RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
								buses.get(0).voltage.getValue());

						RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

						RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

						ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
								tAdmittances);

						double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

						fx0 = new Array2DRowRealMatrix(mismatches);

						Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances,
								indexes, Admittance.createMadmittance(cAdmittances),
								Admittance.createTadmittance(cAdmittances), radmittances, xadmittances);

						RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

						try {

							RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
							double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
									radmittances, xadmittances, DELTA, X0);
							RealMatrix aj = new Array2DRowRealMatrix(arJacob);
							DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
							X1 = X0.add(DELTA);

							if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
								System.out
										.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
								break innerloop;
							}
						} catch (MatrixDimensionMismatchException | DimensionMismatchException
								| NullArgumentException e) {

							e.printStackTrace();
						}

						if (ModifiedNR.sumMatrix(fx0) < 1E-4) {
							boolean flag3 = true;

							if (flag3) {

								flag = false;

							}

						} else if (prev_Mismatches + 1e2 <= ModifiedNR.sumMatrix(fx0)) {

							flag = false;
						} else {

							prev_Mismatches = ModifiedNR.sumMatrix(fx0);

							wi = X1.getEntry(X1.getRowDimension() - 2, 0);

							buses.get(0).voltage = new DerivativeStructure(params, order, 1,
									X1.getEntry(X1.getRowDimension() - 1, 0));

							ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
						}

					}

					double fitness = ModifiedNR.evaluate3(buses, PQLossLoad, wi, position, fx0);

					particles[f].sol.fitness.add(fitness);

				}
			}
			particles[f].setFitness(particles[f].sol.getMean());

		}
	}

	
public Particle[] initializeParticles(int npar) {

		Particle[] particles = new Particle[npar];
		Random random = new Random();
		for (int i = 0; i < npar; i++) {
			double[] initial5 = this.construct2();
			double[] speed = new double[initial5.length];
			for (int j = 0; j < speed.length; j++) {
				speed[j] = random.nextDouble();
			}

			particles[i] = new Particle(initial5, speed);
			particles[i].sol = new Solution2();
			this.experimentParticle(particles[i], 5);

			particles[i].setBestPosition(particles[i].getPosition());
			particles[i].setBestFitness(particles[i].getFitness());

		}

		return particles;
	}

	
public void experimentParticle(Particle p, int rep) {

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;
		RealMatrix prevX1 = null;
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;
		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();

		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadSD);

		p.sol.fitness = new ArrayList<Double>();

		for (int k = 0; k < rep; k++) {
			int params = 12;
			int order = 2;

			X1 = null;
			fx0 = null;
			wi = 1.0;
			prevfx0 = null;
			prevX1 = new Array2DRowRealMatrix(params, 1);
			prev_Mismatches = Double.MAX_VALUE;
			flag = true;
			boolean flag2 = true;
			Jacobian = null;
			double[] position = p.getPosition();

			while (flag2) {
				buses = new ArrayList<Bus>();
				buses.add(new Bus(1, 0, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(1, 2, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(1, 4, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(2, 6, params, 0, 0, position[0], position[2], 0.6, 0.8));
				buses.add(new Bus(1, 8, params, generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 10, params, 0, 0, position[1], position[3], 0.9, 1.2));
				flag2 = ModifiedNR.checkAdequacy(buses, 2.5);
			}
			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 200 && flag; i++) {
				double w0 = 1.00;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

				PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

				RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

				RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

				RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

				ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
						tAdmittances);

				double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

				fx0 = new Array2DRowRealMatrix(mismatches);

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {

					RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
					double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
							radmittances, xadmittances, DELTA, X0);
					RealMatrix aj = new Array2DRowRealMatrix(arJacob);
					DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
					X1 = X0.add(DELTA);

					if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
						System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
						break innerloop;
					}
				} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

					e.printStackTrace();
				}

				if (ModifiedNR.sumMatrix(fx0) < 1E-4) {

					flag = false;

				} else if (prev_Mismatches + 1e2 <= ModifiedNR.sumMatrix(fx0)) {

					X1 = prevX1;
					flag = false;
				} else {

					prevX1 = X1;
					prevfx0 = fx0;
					prev_Mismatches = ModifiedNR.sumMatrix(fx0);

					wi = X1.getEntry(X1.getRowDimension() - 2, 0);

					buses.get(0).voltage = new DerivativeStructure(params, order, 1,
							X1.getEntry(X1.getRowDimension() - 1, 0));

					ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
				}

			}

			// Fitness calculation and update of the fitness ArrayList

			double fitness = ModifiedNR.evaluate3(buses, PQLossLoad, wi, position, fx0);

			p.sol.fitness.add(fitness);

		}

		p.setFitness(p.sol.getMean());

	}

	public double[] construct2() {

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;
		NormalDistribution normal = new NormalDistribution(loadMean, loadSD);
		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;
		boolean flag2 = true;
		double wi = 1.0;
		RealMatrix prevX1 = null;
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;
		Random random = new Random();
		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();

		double[][] Jacobian = null;
		double[] initialPosition = null;
		while (flag2) {
			int params = 12;
			int order = 2;

			X1 = null;
			fx0 = null;
			wi = 1.0;
			prevfx0 = null;
			prevX1 = new Array2DRowRealMatrix(params, 1);
			prev_Mismatches = Double.MAX_VALUE;
			flag = true;
			boolean flag3 = true;

			Jacobian = null;
			while (flag3) {
				buses = new ArrayList<Bus>();
				buses.add(new Bus(1, 0, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(1, 2, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(1, 4, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(2, 6, params, 0, 0, random.nextDouble(), random.nextDouble(), 0.6, 0.8));
				buses.add(new Bus(1, 8, params, generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 10, params, 0, 0, random.nextDouble(), random.nextDouble(), 0.9, 1.2));
				flag3 = ModifiedNR.checkAdequacy(buses, 2.5);

			}

			int nofDroop = 0;
			for (Bus b : buses) {
				if (b.type == 2)
					nofDroop++;
			}

			initialPosition = new double[nofDroop * 2];

			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 200 && flag; i++) {
				double w0 = 1.00;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

				PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

				RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

				RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

				RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

				ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
						tAdmittances);

				double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

				fx0 = new Array2DRowRealMatrix(mismatches);

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {
					try {

						RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
						double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
								radmittances, xadmittances, DELTA, X0);

						RealMatrix aj = new Array2DRowRealMatrix(arJacob);
						DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
						X1 = X0.add(DELTA);

						if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
							System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
							break innerloop;
						}
					} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

						e.printStackTrace();
					}
				} catch (SingularMatrixException e) {

					for (int f = 0; f < Jacobian.length; f++) {
						for (int j = 0; j < Jacobian[1].length; j++) {
							System.out.printf("%.2f \t", Jacobian[f][j]);
						}
						System.out.println();
					}
					for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
						for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
							System.out.printf("%.2f \t", Admittance.createMadmittance(cAdmittances).getEntry(f, j));
						}
						System.out.println();
					}

					for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
						for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
							System.out.printf("%.2f,%.2f \t", cAdmittances[f][j].getReal(),
									cAdmittances[f][j].getImaginary());
						}
						System.out.println();
					}

					break innerloop;

				}
				if (ModifiedNR.sumMatrix(fx0) <= 1E-3) {

//					System.out.println("\nFull convergence achieved. Exiting..." + "*".repeat(10));
					flag = false;
					flag2 = false;

					int y = 0;
					for (Bus b : buses)
						if (b.type == 2)
							initialPosition[y++] = b.mp;

					for (Bus b : buses)
						if (b.type == 2)
							initialPosition[y++] = b.nq;

				} else if (prev_Mismatches + 0.5e2 <= ModifiedNR.sumMatrix(fx0)) {

//					System.out.printf("\nMissed the local optimum. Exiting... \t At iteration : %d", i);
//					System.out.printf("\nRejected sum of mismatches: %.5f", ModifiedNR.sumMatrix(fx0));
					X1 = prevX1;
					flag = false;
				} else {

					prevX1 = X1;
					prevfx0 = fx0;
					prev_Mismatches = ModifiedNR.sumMatrix(fx0);

					wi = X1.getEntry(X1.getRowDimension() - 2, 0);

					buses.get(0).voltage = new DerivativeStructure(params, order, 1,
							X1.getEntry(X1.getRowDimension() - 1, 0));

					ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
				}

			}

		}

		return initialPosition;

	}

	public void addSamples(Particle p, int rep) {

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;
		RealMatrix prevX1 = null;
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();

		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadSD);

		for (int k = 0; k < rep; k++) {
			int params = 12;
			int order = 2;

			X1 = null;
			fx0 = null;
			wi = 1.0;
			prevfx0 = null;
			prevX1 = new Array2DRowRealMatrix(params, 1);
			prev_Mismatches = Double.MAX_VALUE;
			flag = true;
			boolean flag2 = true;
			Jacobian = null;

			double[] position = p.getPosition();

			while (flag2) {
				buses = new ArrayList<Bus>();
				buses.add(new Bus(1, 0, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(1, 2, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(1, 4, params, Math.min(-normal.sample(), 0)));
				buses.add(new Bus(2, 6, params, 0, 0, position[0], position[2], 0.6, 0.8));
				buses.add(new Bus(1, 8, params, generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75)));
				buses.add(new Bus(2, 10, params, 0, 0, position[1], position[3], 0.9, 1.2));
				flag2 = ModifiedNR.checkAdequacy(buses, 2.5);
			}

			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 200 && flag; i++) {
				double w0 = 1.00;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

				PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

				RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

				RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

				RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

				ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
						tAdmittances);

				double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

				fx0 = new Array2DRowRealMatrix(mismatches);

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {

					RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
					double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
							radmittances, xadmittances, DELTA, X0);
					RealMatrix aj = new Array2DRowRealMatrix(arJacob);
					DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
					X1 = X0.add(DELTA);

					if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
						System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
						break innerloop;
					}
				} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

					e.printStackTrace();
				}

				if (ModifiedNR.sumMatrix(fx0) < 1E-4) {
					boolean flag3 = true;

					if (flag3) {

						flag = false;

					}

				} else if (prev_Mismatches + 1e2 <= ModifiedNR.sumMatrix(fx0)) {

					X1 = prevX1;
					flag = false;
				} else {

					prevX1 = X1;
					prevfx0 = fx0;
					prev_Mismatches = ModifiedNR.sumMatrix(fx0);

					wi = X1.getEntry(X1.getRowDimension() - 2, 0);

					buses.get(0).voltage = new DerivativeStructure(params, order, 1,
							X1.getEntry(X1.getRowDimension() - 1, 0));

					ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
				}

			}

			double fitness = ModifiedNR.evaluate3(buses, PQLossLoad, wi, position, fx0);

			p.sol.fitness.add(fitness);

		}

		p.setFitness(p.sol.getMean());

	}

	public static double generateFromWind(double v, double vin, double vout, double vrated, double prated) {

		if (v < vin || v > vout) {
			return 0;
		} else if (v >= vin && v <= vrated) {
			return (prated * (v - vin)) / (vrated - vin);

		} else {
			return prated;
		}

	}

}
