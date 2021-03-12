package nRalgorithm;

public interface FitnessFunction {
	
	public double getFitness(double[] particlePosition);

	public double getFitness(Particle p,int rep);
}
