package nRalgorithm;

public class Main3 {

	public static void main(String[] args) {
		double [] d = new double [] {1,2,3};
		
		double [] d2 = d.clone();
		
		d2[1] = 1000;
		
		System.out.println(d[1]);
		
	}

}
