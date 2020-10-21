package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

public class GenericNR {
	
	public static void main(String[] Args) {
		ArrayList<Double> d = new ArrayList<Double>();
		d.add(1.0);
		d.add(2.0);
		
		ArrayList <String> s = new ArrayList<String>();
		s.add("fp".concat(Integer.toString(1)));

		HashMap<Integer,String> hm = new HashMap<Integer,String>();
		hm.put(1, "fp1");
		System.out.println(hm.get(1));
		double[][] admittances = { { 11, 12, 13 }, { 21, 22, 23 }, { 31, 32, 33 } };
		double[][] theta = { { 11, 12, 13 }, { 21, 22, 23 }, { 31, 32, 33 } };
		int params = 6;
		int order = 2;
		Bus[] buses = new Bus[3];
		buses[0] = new Bus(3, admittances, theta, 0, 0, 0);
		buses[1] = new Bus(0, admittances, theta, 2, -0.9, -0.5);
		buses[2] = new Bus(1, admittances, theta, 4, 0.6, 0);
		buses[0].voltage = new DerivativeStructure(params, order, 1, 1.0);
		buses[0].delta = new DerivativeStructure(params, order, 0, 0);
		buses[2].voltage = new DerivativeStructure(params, order, 5, 1.01);
		System.out.println(buses[0]);
		
		for(int i = 0;i< buses.length;i++) {
			
			
		}
		
		
		
	}
	
	

}
