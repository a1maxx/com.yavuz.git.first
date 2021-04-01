package nRalgorithm;

import java.util.ArrayList;
import java.util.Random;

public class MicrogridFitnessFunction2 implements FitnessFunction {
	@Override
	public double getFitness(double[] particlePosition) {

		double[] mp = new double[3];
		double[] nq = new double[3];
		Random random = new Random();
		
		System.arraycopy(particlePosition, 0, mp, 0, 3);
		System.arraycopy(particlePosition, 3, nq, 0, 3);
		ModifiedNR4 nr4 = new ModifiedNR4();
		ArrayList<Bus> buses = nr4.createBuses(mp, nq);

		return nr4.runMNR(buses);

	}

	@Override
	public double[] getFitness(double[] particlePosition, int rep) {
		double[] mp = new double[3];
		double[] nq = new double[3];
		System.arraycopy(particlePosition, 0, mp, 0, 3);
		System.arraycopy(particlePosition, 3, nq, 0, 3);
		ModifiedNR4 nr4 = new ModifiedNR4();
		ArrayList<Bus> buses = nr4.createBuses(mp, nq);
		double[] sol = new double[rep];
		for (int i = 0; i < rep; i++)
			sol[i] = nr4.runMNR(buses);

		return sol;
	}

	@Override
	public double getFitness(Particle p, int rep) {
		// TODO Auto-generated method stub
		return 0;
	}



}
