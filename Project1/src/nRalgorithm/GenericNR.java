package nRalgorithm;

import java.util.ArrayList;
import java.util.Arrays;


public class GenericNR {
	
	public static void main(String[] Args) {
	

	}
	
	static ArrayList<Integer[]> createOrders(ArrayList<Bus> buses){	
		
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
	return orders;	
	}

}
