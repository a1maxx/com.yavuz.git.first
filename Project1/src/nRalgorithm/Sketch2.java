package nRalgorithm;



//import java.util.Random;

public class Sketch2 {

	public static void main(String[] args) {
		Multiswarm multiswarm = new Multiswarm(5,15,new MicrogridFitnessFunction());

		for (int i = 0; i < 20; i++) {
			multiswarm.mainLoop(i+1);
		}
		System.out.println("Best fitness found: " + multiswarm.getBestFitness() + "\n" +"At position: "+ "["
				+ multiswarm.getBestPosition()[0] + "," + multiswarm.getBestPosition()[1] + "]");
	
	}

}
