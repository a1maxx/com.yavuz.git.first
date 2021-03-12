package nRalgorithm;

import java.util.ArrayList;

import org.apache.commons.math3.distribution.NormalDistribution;

public class ToyFitnessFunction implements FitnessFunction{
	NormalDistribution normal =  new NormalDistribution(0,5.0);
	
	@Override
	public double getFitness(Particle p,int rep) {
	
		double x1 = p.getPosition()[0];
		double x2 = p.getPosition()[1];
		
		p.sol.fitness = new ArrayList<Double>(rep);
		p.sol.rep = rep;
		
	
		for (int i = 0; i < rep; i++) {
			double value = (4 * Math.pow(x1, 2) - 2.1 * Math.pow(x1, 4) + Math.pow(x1, 6) / 3 + x1 * x2
					- 4 * Math.pow(x2, 2) + 4 * Math.pow(x2, 4) + 1.0316) + normal.sample();

//			System.out.printf("\n x1: %.2f \t x2 :%.2f \t value : %.2f",x1,x2,value);
			p.sol.fitness.add(value);
			
		}
				
		return p.sol.getMean();
	}


	@Override
	public double getFitness(double[] particlePosition) {
		// TODO Auto-generated method stub
		return 0;
	}

}
