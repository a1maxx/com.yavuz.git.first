package nRalgorithm;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;

public class ModifiedNR8 {
	public static final double loadMean = 0.1;

	double[][] xadmittances;
	double[][] radmittances;
	double[][] theta;

	public ModifiedNR8() {
		int nofB = 6;
		double[][] xadmittances = new double[nofB][nofB];
		double[][] radmittances = new double[nofB][nofB];
		File file = new File("test4.txt");
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

	public static void main(String[] args) {
		ModifiedNR8 mnr8 = new ModifiedNR8();

		int npar = 5;

		Particle[] particles = mnr8.initializeParticles(npar);

		int iteration = 1;

		int max_iter = 40;

		double bestFitness = Double.MAX_VALUE;
		double[] bestPosition = new double[6];

		int n0 = 20;

		bestFitness = particles[ModifiedNR7.findBestInd(particles)].getFitness();
		bestPosition = particles[ModifiedNR7.findBestInd(particles)].getPosition().clone();

		particles = mnr8.initializeParticles(npar);

		for (iteration = 1; iteration < max_iter; iteration++) {
			Particle best = new Particle();
			best.setPosition(bestPosition);
			best.setFitness(bestFitness);
			mnr8.addSamples(best, n0);

			bestFitness = best.getFitness();

			boolean flag = false;
			for (Particle particle : particles) {

				for (Particle particle2 : particles)
					mnr8.experimentParticle(particle2, n0);

//Case1
//				int budget =50;
//				mnr8.makeAdditionalReps0(particles,budget);

//Case 2

//				int budget = 50;
//				mnr8.makeAdditionalReps2(particles, budget, bestFitness);

//////Case 3
				int budget = 50;
				Res res = ModifiedNR7.findCase(particles, bestFitness);
				mnr8.makeAdditionalReps3(particles, budget, bestFitness, res, best);

				double[] pOldPos = particle.getPosition().clone();

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

				// Update the particle speed
				for (int m = 0; m < speed.length; m++) {
					speed[m] = ModifiedNR7.getNewParticleSpeedForIndex(particle, m, iteration, max_iter, bestPosition);
				}

			}
			// If best changes, print
			if (flag) {
				ModifiedNR7.printBest(bestPosition, iteration);
				System.out.printf("%.4f\n", bestFitness);
			} else {
				if (iteration % 100 != 0)
					System.out.print(".");
				else
					System.out.println();
			}
		}

	}

	public void makeAdditionalReps2(Particle[] particles, int budget, double bestFitness) {
		int[] adReps = ModifiedNR7.findRep(particles, budget);

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;
		RealMatrix prevX1 = new Array2DRowRealMatrix(12, 1);
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.1);

		for (int f = 0; f < particles.length; f++) {
			int rep = adReps[f];

			if (rep > 0) {
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

					Jacobian = null;
					double PF = 0.628;

					double[] position = particles[f].getPosition();

					buses = new ArrayList<Bus>();
					buses.add(new Bus(1, 0, params, -normal.sample()));
					buses.add(new Bus(1, 2, params, 0.05, 0.05 * PF));
					buses.add(new Bus(1, 8, params, -normal.sample()));
					buses.add(new Bus(2, 4, params, 0, 0, position[0], position[3], 3.0));
					buses.add(new Bus(2, 6, params, 0, 0, position[1], position[4], 3.0));
					buses.add(new Bus(2, 10, params, 0, 0, position[2], position[5], 3.0));

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

					double fitness = Math.abs(PQLossLoad.get(0)) + Math.abs(PQLossLoad.get(1));
					for (Bus b : buses)
						fitness = fitness + Math.abs(1 - b.voltage.getValue());

					if (ModifiedNR.sumMatrix(fx0) > 1E-4)
						fitness += ModifiedNR.sumMatrix(fx0);

					fitness += Math.abs(1 - wi);

					particles[f].sol.fitness.add(fitness);

				}
			}
			particles[f].setFitness(particles[f].sol.getMean());

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

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.1);

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
					double PF = 0.628;

					double[] position = particles[f].getPosition();

					buses = new ArrayList<Bus>();
					buses.add(new Bus(1, 0, params, -normal.sample()));
					buses.add(new Bus(1, 2, params, 0.05, 0.05 * PF));
					buses.add(new Bus(1, 4, params, -normal.sample()));
					buses.add(new Bus(2, 6, params, 0, 0, position[0], position[3], 3.0));
					buses.add(new Bus(2, 8, params, 0, 0, position[1], position[4], 3.0));
					buses.add(new Bus(2, 10, params, 0, 0, position[2], position[5], 3.0));

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

					double fitness = Math.abs(PQLossLoad.get(0)) + Math.abs(PQLossLoad.get(1));
					for (Bus b : buses) {
						fitness = fitness + Math.abs(1 - b.voltage.getValue());

					}

					if (ModifiedNR.sumMatrix(fx0) > 1E-4)
						fitness += ModifiedNR.sumMatrix(fx0);

					fitness += Math.abs(1 - wi);

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
			double[] initial5 = this.construct();
			double[] speed = new double[initial5.length];
			for (int j = 0; j < speed.length; j++) {
				speed[j] = random.nextDouble();
			}

			particles[i] = new Particle(initial5, speed);
			particles[i].sol = new Solution2();
			this.experimentParticle(particles[i], 10);

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

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.1);

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

			Jacobian = null;
			double PF = 0.628;

			double[] position = p.getPosition();
			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, params, -normal.sample()));
			buses.add(new Bus(1, 2, params, 0.05, 0.05 * PF));
			buses.add(new Bus(1, 4, params, -normal.sample()));
			buses.add(new Bus(2, 6, params, 0, 0, position[0], position[3], 3.0));
			buses.add(new Bus(2, 8, params, 0, 0, position[1], position[4], 3.0));
			buses.add(new Bus(2, 10, params, 0, 0, position[2], position[5], 3.0));

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

			double fitness = Math.abs(PQLossLoad.get(0)) + Math.abs(PQLossLoad.get(1));
			for (Bus b : buses) {
				fitness = fitness + Math.abs(1 - b.voltage.getValue());

			}

			if (ModifiedNR.sumMatrix(fx0) > 1E-4)
				fitness += ModifiedNR.sumMatrix(fx0);

			fitness += Math.abs(1 - wi);

			p.sol.fitness.add(fitness);

		}

		p.setFitness(p.sol.getMean());

	}

	public double[] construct() {

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;
		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.1);
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

			Jacobian = null;
			double PF = 0.628;

			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, params, -normal.sample()));
			buses.add(new Bus(1, 2, params, 0.05, 0.05 * PF));
			buses.add(new Bus(1, 4, params, -normal.sample()));
			buses.add(new Bus(2, 6, params, 0, 0, (random.nextDouble() * 2) - 1, (random.nextDouble() * 2), 3.0));
			buses.add(new Bus(2, 8, params, 0, 0, (random.nextDouble() * 2) - 1, (random.nextDouble() * 2), 3.0));
			buses.add(new Bus(2, 10, params, 0, 0, (random.nextDouble() * 2) - 1, (random.nextDouble() * 2), 3.0));

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
				if (ModifiedNR.sumMatrix(fx0) <= 1E-4) {

					System.out.println("\nFull convergence achieved. Exiting..." + "*".repeat(10));
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

					System.out.printf("\nMissed the local optimum. Exiting... \t At iteration : %d", i);
					System.out.printf("\nRejected sum of mismatches: %.5f", ModifiedNR.sumMatrix(fx0));
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

	public void makeAdditionalReps0(Particle[] particles, int budget) {
		int[] adReps = new int[particles.length];
		Arrays.fill(adReps, budget / particles.length);

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;
		RealMatrix prevX1 = new Array2DRowRealMatrix(12, 1);
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.1);

		for (int f = 0; f < particles.length; f++) {
			int rep = adReps[f];

			if (rep > 0) {
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

					Jacobian = null;
					double PF = 0.628;

					double[] position = particles[f].getPosition();

					buses = new ArrayList<Bus>();
					buses.add(new Bus(1, 0, params, -normal.sample()));
					buses.add(new Bus(1, 2, params, 0.05, 0.05 * PF));
					buses.add(new Bus(1, 4, params, -normal.sample()));
					buses.add(new Bus(2, 6, params, 0, 0, position[0], position[3], 3.0));
					buses.add(new Bus(2, 8, params, 0, 0, position[1], position[4], 3.0));
					buses.add(new Bus(2, 10, params, 0, 0, position[2], position[5], 3.0));

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

					double fitness = Math.abs(PQLossLoad.get(0)) + Math.abs(PQLossLoad.get(1));

					for (Bus b : buses)
						fitness = fitness + Math.abs(1 - b.voltage.getValue());

					if (ModifiedNR.sumMatrix(fx0) > 1E-4)
						fitness += ModifiedNR.sumMatrix(fx0);

					fitness += Math.abs(1 - wi);

					particles[f].sol.fitness.add(fitness);

				}
			}
			particles[f].setFitness(particles[f].sol.getMean());

		}
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
		Random random = new Random();
		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		long cur = System.currentTimeMillis();
		double[][] Jacobian = null;

		NormalDistribution normal = new NormalDistribution(loadMean, loadMean * 0.1);

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

			Jacobian = null;
			double PF = 0.628;

			double[] position = p.getPosition();
			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, params, -normal.sample()));
			buses.add(new Bus(1, 2, params, 0.05, 0.05 * PF));
			buses.add(new Bus(1, 4, params, -normal.sample()));
			buses.add(new Bus(2, 6, params, 0, 0, position[0], position[3], 3.0));
			buses.add(new Bus(2, 8, params, 0, 0, position[1], position[4], 3.0));
			buses.add(new Bus(2, 10, params, 0, 0, position[2], position[5], 3.0));

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

			double fitness = Math.abs(PQLossLoad.get(0)) + Math.abs(PQLossLoad.get(1));
			for (Bus b : buses) {
				fitness = fitness + Math.abs(1 - b.voltage.getValue());

			}

			if (ModifiedNR.sumMatrix(fx0) > 1E-4)
				fitness += ModifiedNR.sumMatrix(fx0);

			fitness += Math.abs(1 - wi);

			p.sol.fitness.add(fitness);

		}

		p.setFitness(p.sol.getMean());

	}
}
