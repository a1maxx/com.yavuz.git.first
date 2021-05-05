package nRalgorithm;
import java.io.File;
import java.io.IOException;

import java.text.SimpleDateFormat;
import java.util.Date;

public class Main3 {

	public static void main(String[] args) throws IOException  {
		
		
	      try {
	    	  String fileName = new SimpleDateFormat("yyyyMMddHHmmss'.txt'").format(new Date());
	    	  fileName = "zoutput_" + fileName ;
	          File myObj = new File(fileName);
	          if (myObj.createNewFile()) {
	            System.out.println("File created: " + myObj.getName());
	          } else {
	            System.out.println("File already exists.");
	          }
	        } catch (IOException e) {
	          System.out.println("An error occurred.");
	          e.printStackTrace();
	        }


		
	}

}
