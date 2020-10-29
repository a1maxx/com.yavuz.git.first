package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;

public class GenericNR {
	
	public static void main(String[] Args) {
		ArrayList<Bus> buses = new ArrayList<Bus>();

		ArrayList <String> s = new ArrayList<String>();
		s.add("fp".concat(Integer.toString(1)));
		
		double[][] admittances = { { 11, 12, 13 }, { 21, 22, 23 }, { 31, 32, 33 } };
		double[][] theta = { { 11, 12, 13 }, { 21, 22, 23 }, { 31, 32, 33 } };
		int params = 6;
		int order = 2;
		
		buses.add(new Bus(3, admittances, theta, 0, 0, 0));
		buses.add(new Bus(0, admittances, theta, 2, -0.9, -0.5));
		buses.add(new Bus(1, admittances, theta, 4, 0.6, 0));
		buses.get(0).voltage = new DerivativeStructure(params, order, 1, 1.0);
		buses.get(0).delta = new DerivativeStructure(params, order, 0, 0);
		buses.get(2).voltage = new DerivativeStructure(params, order, 5, 1.01);
		//System.out.println(buses);
		
		ArrayList<Integer[]> orders = new ArrayList<Integer []>();
		
		for(Bus b : buses ) {
			
			switch (b.type) {
				case 0:
					Integer [] temp = new Integer[buses.size()*2];
					Arrays.fill(temp,0);
					temp[b.index] = 1;
					orders.add(temp);
					temp = new Integer[buses.size()*2];
					Arrays.fill(temp, 0);
					temp[b.index+1] = 1;
					orders.add(temp);
					break;
				case 1:
					temp = new Integer[buses.size()*2];
					Arrays.fill(temp,0);
					temp[b.index] = 1;
					orders.add(temp);
					break;
				case 2:
					temp = new Integer[buses.size()*2];
					Arrays.fill(temp, 0);
					temp[b.index+1] = 1;
					orders.add(temp);
					break;
				default:

			}
			
		}
		Bus.createEqauations2(buses);
	

	}
	
	

}
