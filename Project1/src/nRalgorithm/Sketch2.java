package nRalgorithm;

import org.apache.commons.math3.distribution.WeibullDistribution;

import java.util.Random;

public class Sketch2 {

	public static void main(String[] args) {
		WeibullDistribution wb2 = new WeibullDistribution(7.5, 3.5);
		WeibullDistribution wb = new WeibullDistribution(7.5, 7);
		Random r =  new Random();
		for(int i=0; i< 100; i++) {
			System.out.printf("WB1:\t %.3f \tWB2:\t %.3f \n",ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75),
			ModifiedNR8.generateFromWind(wb2.sample(), 3.5, 20.0, 14.5, 0.75));

		}

		
//		Multiswarm multiswarm = new Multiswarm(10,25 ,new MicrogridFitnessFunction2());
//
//		for (int i = 0; i < 200; i++) {
//			multiswarm.mainLoop(i+1);
//		}
//		
//		
//		
//		System.out.println("Best fitness found: " + multiswarm.getBestFitness());
//		
//		for(int i=0; i < multiswarm.getBestPosition().length;i++) {
//			System.out.printf("%.8f, ",multiswarm.getBestPosition()[i]);
//			
//		}
//		System.out.println();
		
		
	}

}
