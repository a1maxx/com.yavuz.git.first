package nRalgorithm;



//import java.util.Random;

public class Sketch2 {

	public static void main(String[] args) {
		Multiswarm multiswarm = new Multiswarm(10,25 ,new MicrogridFitnessFunction2());

		for (int i = 0; i < 200; i++) {
			multiswarm.mainLoop(i+1);
		}
		
		
		
		System.out.println("Best fitness found: " + multiswarm.getBestFitness());
		
		for(int i=0; i < multiswarm.getBestPosition().length;i++) {
			System.out.printf("%.8f, ",multiswarm.getBestPosition()[i]);
			
		}
		System.out.println();
		
		
	}

}
