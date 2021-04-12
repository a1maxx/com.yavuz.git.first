package nRalgorithm;


import java.util.Random;

import org.apache.commons.math3.distribution.WeibullDistribution;

public class Main2 {

	public static void main(String[] args) {

		Random random = new Random();
		
//		Test0 tst = new Test0();
//		tst.printTo(null, null);
//		
		WeibullDistribution wb = new WeibullDistribution(7.5, 3.25);
		double v = 0;
		double generation = 0;
		int count = 0;
		for (int i = 0; i < 100; i++) {
			v = wb.sample();
//			System.out.println(v);
			generation = ModifiedNR0.generateFromWind(v, 3.5, 20.0, 14.5, 0.75);
			if (generation > 0) {
				count++;
				System.out.println(generation);
			} 
			
		
		}
		System.out.println(count);
//		
//		BetaDistribution bt = new BetaDistribution(0.40,8.56);
//		for(int i =0; i<100;i++)
//			System.out.println(bt.sample());
		

		
//		int nPar = 10;
//		Particle particles[] = new Particle[nPar];
//		for (int i = 0; i < nPar; i++) {
//			particles[i] = new Particle(
//					new double[] { 200 * random.nextDouble() - 100, 200 * random.nextDouble() - 100 },
//					new double[] { 200 * random.nextDouble() - 100, 200 * random.nextDouble() - 100 });
//
//			particles[i].sol = new Solution2();
//			
//		}
//		ToyFitnessFunction tff = new ToyFitnessFunction();	
//		for(Particle p: particles) {
//			p.setFitness(tff.getFitness(p, p.getRep()));
//			if(p.getFitness()< p.getBestFitness())
//				p.setBestFitness(p.getFitness());
//			
//		}
//		double currentBest = particles[ModifiedNR7.findBestInd(particles)].getFitness();
//		
//		for(Particle p: particles)
//			p.setFitness(tff.getFitness(p, p.getRep()));
//		

//		ModifiedNR7.findRep2(particles, 100,currentBest);
//		
		
		
//		DerivativeStructure d1 = new DerivativeStructure(2,2,0,10) ;
//		DerivativeStructure d2 = new DerivativeStructure(2,2,1,20) ;
//	
//	
//		DerivativeStructure d3 =  new DerivativeStructure(2,2);
//			
//		d3 = d1.pow(2).add(d2) ;
//		
//		System.out.println(d3.getPartialDerivative(1,0));
//		
//		d1 = new DerivativeStructure(2,2,0,5) ;
//		
//		d3 = d1.pow(2).add(d2) ;
//		
//		System.out.println(d3.getPartialDerivative(1,0));

	}
	


}
