package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

public class GenericNR {
	
	public static void main(String[] Args) {
		ArrayList<Double> d = new ArrayList<Double>();
		d.add(1.0);
		d.add(2.0);
		System.out.println(d.get(0));
		System.out.println(d.indexOf(1.0));
		ArrayList <String> s = new ArrayList<String>();
		s.add("fp".concat(Integer.toString(1)));
		System.out.print(s.indexOf("fp".concat(Integer.toString(1))));
		System.out.print("\na");
		
	}
	
	

}
