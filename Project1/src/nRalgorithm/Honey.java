package nRalgorithm;
/* Honey.java
 *
 * Honey class used by ArtificialBeeColony.java
 * Also known as food source.
 * Contains the positions of the queens in a solution as well as its conflicts, trials, fitness and selection probability. 
 * Base code at http://mf.erciyes.edu.tr/abc/
 *
 * @author: James M. Bayon-on
 * @version: 1.3
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.lang.Double;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;



public class Honey implements Comparable<Honey> {
	private int MAX_LENGTH;
	private double nectar[]; 							//solution or placement of queens
    private int trials;
    private double conflicts;
    private double fitness;
    private double selectionProbability;
	double[][] xadmittances;
	double[][] radmittances;
	double[][] theta;
	public static final double SENSITIVITY = 1E-3;
	
	WeibullDistribution wb = new WeibullDistribution(7.5, 3.5);
	double 	wGen = ModifiedNR8.generateFromWind(wb.sample(), 3.5, 20.0, 14.5, 0.75);
    /* Instantiate a Honey.
     *
     * @param: size of n
     */
    public Honey(int size) {
        this.MAX_LENGTH = size;
        this.nectar = new double[MAX_LENGTH];
        this.conflicts = 0;
        this.trials = 0;
        this.fitness = 0.0;
        this.selectionProbability = 0.0;
//		int nofB = 33;
//		double[][] xadmittances = new double[nofB][nofB];
//		double[][] radmittances = new double[nofB][nofB];
//		File file = new File("case33.txt");
//		try {
//			Scanner scan = new Scanner(file);
//			for (int i = 0; i < xadmittances.length; i++) {
//				for (int j = 0; j < xadmittances[1].length; j++) {
//					xadmittances[i][j] = scan.nextDouble();
//
//				}
//			}
//			for (int i = 0; i < xadmittances.length; i++) {
//				for (int j = 0; j < xadmittances[1].length; j++) {
//					radmittances[i][j] = scan.nextDouble();
//				}
//			}
//			scan.close();
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
		int nofB = 6;
		double[][] xadmittances = new double[nofB][nofB];
		double[][] radmittances = new double[nofB][nofB];
		File file = new File("case6ww.txt");
		try {
			Scanner scan = new Scanner(file);
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					xadmittances[i][j] = scan.nextDouble();

				}
			}
			for (int i = 0; i < xadmittances.length; i++) {
				for (int j = 0; j < xadmittances[1].length; j++) {
					radmittances[i][j] = scan.nextDouble();
				}
			}
			scan.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		this.xadmittances = xadmittances;
		this.radmittances = radmittances;

		this.xadmittances = xadmittances;
		this.radmittances = radmittances;
        initNectar();
    }
    
    /* Compares two Honeys.
	 *
	 * @param: a Honey to compare with
	 */	
//    public int compareTo0(Honey h) {
//   
//    	return this.conflicts - h.getConflicts();
//    }
    
  public int compareTo(Honey h) {
  	return Double.compare(this.conflicts, h.getConflicts());
 }

    /* Initializes the Honey into diagonal queens.
	 *
	 */
    public void initNectar0() { 
        for(int i = 0; i < MAX_LENGTH; i++) {		//initialize the solution to 1... n
            nectar[i] = i;
        }
    }
    
    public void initNectar() {
    	double [] initial = this.construct();
        for(int i = 0; i < MAX_LENGTH; i++) {		//initialize the solution to 1... n
            nectar[i] = initial[i];
        }
    }
    
    
    
    public void computeConflicts0() { 				//compute the number of conflicts to calculate fitness

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;

		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;
		
		X1 = null;
		fx0 = null;
		wi = 1.0;


		int params = 66;
		int order = 2;

		Jacobian = null;

		double[] position = this.nectar;

		buses = new ArrayList<Bus>();
		buses = new ArrayList<Bus>();
		buses.add(new Bus(1, 0, params,  Math.min(-0.0015,0)));//1
		buses.add(new Bus(1, 2, params,  Math.min(-0.0015,0)));//2
		buses.add(new Bus(1, 4, params,  Math.min(-0.0015,0)));//3
		buses.add(new Bus(1, 6, params,  Math.min(-0.0015,0)));//4
		buses.add(new Bus(1, 8, params,  Math.min(-0.0015,0)));//5
		buses.add(new Bus(1, 10, params, wGen)); // 6
		buses.add(new Bus(1, 12, params,  Math.min(-0.0015,0)));//7
		buses.add(new Bus(1, 14, params,  Math.min(-0.0015,0)));//8
		buses.add(new Bus(1, 16, params,  Math.min(-0.0015,0)));//9
		buses.add(new Bus(1, 18, params, wGen)); //10
		buses.add(new Bus(1, 20, params,  Math.min(-0.0015,0)));//11
		buses.add(new Bus(1, 22, params,  Math.min(-0.0015,0)));//12
		buses.add(new Bus(2, 24, params, 0, 0, position[0], position[4], 1.8, 2.4)); //13
		buses.add(new Bus(1, 26, params,  Math.min(-0.0015,0)));//14
		buses.add(new Bus(2, 28, params, 0, 0, position[1], position[5], 1.2, 1.6)); //15
		buses.add(new Bus(1, 30, params,  Math.min(-0.0015,0)));//16
		buses.add(new Bus(1, 32, params,  Math.min(-0.0015,0)));//17
		buses.add(new Bus(1, 34, params,  Math.min(-0.0015,0)));//18
		buses.add(new Bus(1, 36, params,  Math.min(-0.0015,0)));//19
		buses.add(new Bus(1, 38, params,  wGen));//20
		buses.add(new Bus(1, 40, params,  Math.min(-0.0015,0)));//21
		buses.add(new Bus(1, 42, params,  Math.min(-0.0015,0)));//22
		buses.add(new Bus(1, 44, params,  Math.min(-0.0015,0)));//23
		buses.add(new Bus(1, 46, params,  Math.min(-0.0015,0)));//24
		buses.add(new Bus(2, 48, params, 0, 0, position[2], position[6], 1.2, 1.6)); //25
		buses.add(new Bus(1, 50, params,  Math.min(-0.0015,0)));//26
		buses.add(new Bus(1, 52, params,  Math.min(-0.0015,0)));//27
		buses.add(new Bus(1, 54, params,  Math.min(-0.0015,0)));//28
		buses.add(new Bus(1, 56, params,  wGen)); //29
		buses.add(new Bus(1, 58, params,  Math.min(-0.0015,0)));//30
		buses.add(new Bus(1, 60, params,  Math.min(-0.0015,0)));//31
		buses.add(new Bus(1, 62, params,  Math.min(-0.0015,0)));//32
		buses.add(new Bus(2, 64, params, 0, 0, position[3], position[7], 1.2, 1.6)); //33

		
		
		
		
		
		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

		innerloop: for (int i = 1; i <= 200 && flag; i++) {
			double w0 = 1.00;
			double v0 = 1.01;

			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances,
					xadmittances, wi);

			ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

			PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

			RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
					buses.get(0).voltage.getValue());

			RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

			RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
					tAdmittances);

			double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

			fx0 = new Array2DRowRealMatrix(mismatches);

			Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances,
					indexes, Admittance.createMadmittance(cAdmittances),
					Admittance.createTadmittance(cAdmittances), radmittances, xadmittances);

			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

			try {

				RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
				double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
						radmittances, xadmittances, DELTA, X0);
				RealMatrix aj = new Array2DRowRealMatrix(arJacob);
				DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
				X1 = X0.add(DELTA);

				if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
					System.out
							.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
					break innerloop;
				}
			} catch (MatrixDimensionMismatchException | DimensionMismatchException
					| NullArgumentException e) {

				e.printStackTrace();
			}

			if (ModifiedNR.sumMatrix(fx0) < 1E-4) {
				boolean flag3 = true;

				if (flag3) {

					flag = false;

				}

			} else if (prev_Mismatches + 1e2 <= ModifiedNR.sumMatrix(fx0)) {

				flag = false;
			} else {

				prev_Mismatches = ModifiedNR.sumMatrix(fx0);

				wi = X1.getEntry(X1.getRowDimension() - 2, 0);

				buses.get(0).voltage = new DerivativeStructure(params, order, 1,
						X1.getEntry(X1.getRowDimension() - 1, 0));

				ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
			}

		}

		
		
		this.conflicts= ModifiedNR.evaluate3(buses, PQLossLoad, wi, position, fx0);;
    }
    
    
    
    public void computeConflicts() { 				//compute the number of conflicts to calculate fitness
    	int params = 12;
		int order = 2;
		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;

		double wi = 1.0;

		double prev_Mismatches = Double.MAX_VALUE;

		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();
		double[][] Jacobian = null;
		
		X1 = null;
		fx0 = null;
		wi = 1.0;

		prev_Mismatches = Double.MAX_VALUE;
		flag = true;

		Jacobian = null;

		double[] position = this.nectar;

		buses = new ArrayList<Bus>();
		buses.add(new Bus(1, 0, params, -0.1));
		buses.add(new Bus(1, 2, params, 0.01));
		buses.add(new Bus(1, 4, params, -0.1));
		buses.add(new Bus(2, 6, params, 0, 0, position[0], position[3], 0.9,1.2));
		buses.add(new Bus(2, 8, params, 0, 0, position[1], position[4], 0.9,1.2));
		buses.add(new Bus(2, 10, params, 0, 0, position[2], position[5], 0.9,1.2));

		ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
		ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

		innerloop: for (int i = 1; i <= 200 && flag; i++) {
			double w0 = 1.00;
			double v0 = 1.01;

			Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances,
					xadmittances, wi);

			ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

			PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

			RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
					buses.get(0).voltage.getValue());

			RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

			RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

			ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
					tAdmittances);

			double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

			fx0 = new Array2DRowRealMatrix(mismatches);

			Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances,
					indexes, Admittance.createMadmittance(cAdmittances),
					Admittance.createTadmittance(cAdmittances), radmittances, xadmittances);

			RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

			try {

				RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
				double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
						radmittances, xadmittances, DELTA, X0);
				RealMatrix aj = new Array2DRowRealMatrix(arJacob);
				DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
				X1 = X0.add(DELTA);

				if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
					System.out
							.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
					break innerloop;
				}
			} catch (MatrixDimensionMismatchException | DimensionMismatchException
					| NullArgumentException e) {

				e.printStackTrace();
			}

			if (ModifiedNR.sumMatrix(fx0) < 1E-4) {
				boolean flag3 = true;

				if (flag3) {

					flag = false;

				}

			} else if (prev_Mismatches + 1e2 <= ModifiedNR.sumMatrix(fx0)) {

				flag = false;
			} else {

				prev_Mismatches = ModifiedNR.sumMatrix(fx0);

				wi = X1.getEntry(X1.getRowDimension() - 2, 0);

				buses.get(0).voltage = new DerivativeStructure(params, order, 1,
						X1.getEntry(X1.getRowDimension() - 1, 0));

				ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
			}

		}

		
//		
//		
//		double fitness = Math.abs(PQLossLoad.get(0)) + Math.abs(PQLossLoad.get(1));
//		for (Bus b : buses) {
//			fitness = fitness + Math.abs(1 - b.voltage.getValue());
//
//		}
//
//		if (ModifiedNR.sumMatrix(fx0) > 1E-4)
//			fitness += ModifiedNR.sumMatrix(fx0);
//
//		fitness += (Math.abs(1 - wi) * 10) * buses.size();
//		
		
		this.conflicts= ModifiedNR.evaluate3(buses, PQLossLoad, wi, position, fx0);;
    }
    
    
    
    
    /* Computes the conflicts in the nxn board.
	 *
	 */
//    public void computeConflicts0() { 				//compute the number of conflicts to calculate fitness
//        String board[][] = new String[MAX_LENGTH][MAX_LENGTH]; 
//        int x = 0; 									//row
//        int y = 0; 									//column
//        int tempx = 0;
//        int tempy = 0; 
//        
//        int dx[] = new int[] {-1, 1, -1, 1}; 		//to check for diagonal
//        int dy[] = new int[] {-1, 1, 1, -1}; 		//paired with dx to check for diagonal
//        
//        boolean done = false; 						//used to check is checking fo diagonal is out of bounds
//        int conflicts = 0; 							//number of conflicts found
//        
//        board = clearBoard(board); 					//clears the board into empty strings
//        board = plotQueens(board); 					//plots the Q in the board
// 
//        // Walk through each of the Queens and compute the number of conflicts.
//        for(int i = 0; i < MAX_LENGTH; i++) {
//            x = i;
//            y = this.nectar[i]; 					//will result to no horizontal and vertical conflicts because it will result to diagonal 
//
//            // Check diagonals.
//            for(int j = 0; j < 4; j++) { 			//because of dx and dy where there are 4 directions for diagonal searching for conflicts
//                tempx = x;
//                tempy = y; 							// store coordinate in temp
//                done = false;
//                
//                while(!done) {
//                    tempx += dx[j];
//                    tempy += dy[j];
//                    
//                    if((tempx < 0 || tempx >= MAX_LENGTH) || (tempy < 0 || tempy >= MAX_LENGTH)) { //if exceeds board
//                        done = true;
//                    } else {
//                        if(board[tempx][tempy].equals("Q")) {
//                            conflicts++;
//                        }
//                    }
//                }
//            }
//        }
//
//        this.conflicts = conflicts; 				//set conflicts of this chromosome
//        
//    }
    
    /* Plots the queens in the board.
	 *
	 * @param: a nxn board
	 */
//    public String[][] plotQueens(String[][] board) {
//        for(int i = 0; i < MAX_LENGTH; i++) {
//            board[i][this.nectar[i]] = "Q";
//        }
//        return board;
//    }
    
    /* Clears the board.
	 *
	 * @param: a nxn board
	 */
    public String[][] clearBoard(String[][] board) {
        // Clear the board.
        for(int i = 0; i < MAX_LENGTH; i++) {
            for(int j = 0; j < MAX_LENGTH; j++) {
                board[i][j] = "";
            }
        }
        return board;
    }

    /* Gets the conflicts of the Honey.
	 *
	 * @return: number of conflicts of the honey
	 */
    public double getConflicts() {
        return conflicts;
    }

    
	public double[] construct0() {

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;
		boolean flag2 = true;
		double wi = 1.0;
		RealMatrix prevX1 = null;
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;
		Random random = new Random();
		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();

		double[][] Jacobian = null;
		double[] initialPosition = null;
		while (flag2) {
			int params = 66;
			int order = 2;

			X1 = null;
			fx0 = null;
			wi = 1.0;
			prevfx0 = null;
			prevX1 = new Array2DRowRealMatrix(params, 1);
			prev_Mismatches = Double.MAX_VALUE;
			flag = true;
			boolean flag3 = true;
			Jacobian = null;
			int FAC = 4;
			int MIN = 2;
			
			while (flag3) {
				
				buses = new ArrayList<Bus>();
				buses.add(new Bus(1, 0, params,  Math.min(-0.0015,0)));//1
				buses.add(new Bus(1, 2, params,  Math.min(-0.0015,0)));//2
				buses.add(new Bus(1, 4, params,  Math.min(-0.0015,0)));//3
				buses.add(new Bus(1, 6, params,  Math.min(-0.0015,0)));//4
				buses.add(new Bus(1, 8, params,  Math.min(-0.0015,0)));//5
				buses.add(new Bus(1, 10, params, wGen)); // 6
				buses.add(new Bus(1, 12, params,  Math.min(-0.0015,0)));//7
				buses.add(new Bus(1, 14, params,  Math.min(-0.0015,0)));//8
				buses.add(new Bus(1, 16, params,  Math.min(-0.0015,0)));//9
				buses.add(new Bus(1, 18, params, wGen)); //10
				buses.add(new Bus(1, 20, params,  Math.min(-0.0015,0)));//11
				buses.add(new Bus(1, 22, params,  Math.min(-0.0015,0)));//12
				buses.add(new Bus(2, 24, params, 0, 0, random.nextDouble()*FAC-MIN, random.nextDouble()*FAC-MIN, 1.8, 2.4)); //13
				buses.add(new Bus(1, 26, params,  Math.min(-0.0015,0)));//14
				buses.add(new Bus(2, 28, params, 0, 0, random.nextDouble()*FAC-MIN, random.nextDouble()*FAC-MIN, 1.2, 1.6)); //15
				buses.add(new Bus(1, 30, params,  Math.min(-0.0015,0)));//16
				buses.add(new Bus(1, 32, params,  Math.min(-0.0015,0)));//17
				buses.add(new Bus(1, 34, params,  Math.min(-0.0015,0)));//18
				buses.add(new Bus(1, 36, params,  Math.min(-0.0015,0)));//19
				buses.add(new Bus(1, 38, params,  wGen));//20
				buses.add(new Bus(1, 40, params,  Math.min(-0.0015,0)));//21
				buses.add(new Bus(1, 42, params,  Math.min(-0.0015,0)));//22
				buses.add(new Bus(1, 44, params,  Math.min(-0.0015,0)));//23
				buses.add(new Bus(1, 46, params,  Math.min(-0.0015,0)));//24
				buses.add(new Bus(2, 48, params, 0, 0, random.nextDouble()*FAC-MIN, random.nextDouble()*FAC-MIN, 1.2, 1.6)); //25
				buses.add(new Bus(1, 50, params,  Math.min(-0.0015,0)));//26
				buses.add(new Bus(1, 52, params,  Math.min(-0.0015,0)));//27
				buses.add(new Bus(1, 54, params,  Math.min(-0.0015,0)));//28
				buses.add(new Bus(1, 56, params, wGen)); //29
				buses.add(new Bus(1, 58, params,  Math.min(-0.0015,0)));//30
				buses.add(new Bus(1, 60, params,  Math.min(-0.0015,0)));//31
				buses.add(new Bus(1, 62, params,  Math.min(-0.0015,0)));//32
				buses.add(new Bus(2, 64, params, 0, 0, random.nextDouble()*FAC-MIN, random.nextDouble()*FAC-MIN, 1.2, 1.6)); //33

	
				
				flag3 = ModifiedNR.checkAdequacy(buses, 7);
			}

			int nofDroop = 0;
			for (Bus b : buses) {
				if (b.type == 2)
					nofDroop++;
			}

			initialPosition = new double[nofDroop * 2];

			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 100 && flag; i++) {  
				double w0 = 1.00;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

				PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

				RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

				RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

				RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

				ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
						tAdmittances);

				double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

				fx0 = new Array2DRowRealMatrix(mismatches);

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {
					try {

						RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
						double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
								radmittances, xadmittances, DELTA, X0);

						RealMatrix aj = new Array2DRowRealMatrix(arJacob);
						DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
						X1 = X0.add(DELTA);

						if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
							System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
							break innerloop;
						}
					} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

						e.printStackTrace();
					}
				} catch (SingularMatrixException e) {

					for (int f = 0; f < Jacobian.length; f++) {
						for (int j = 0; j < Jacobian[1].length; j++) {
							System.out.printf("%.2f \t", Jacobian[f][j]);
						}
						System.out.println();
					}
					for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
						for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
							System.out.printf("%.2f \t", Admittance.createMadmittance(cAdmittances).getEntry(f, j));
						}
						System.out.println();
					}

					for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
						for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
							System.out.printf("%.2f,%.2f \t", cAdmittances[f][j].getReal(),
									cAdmittances[f][j].getImaginary());
						}
						System.out.println();
					}

					break innerloop;

				}
				if (ModifiedNR.sumMatrix(fx0) <= SENSITIVITY) {

					System.out.println("\nFull convergence achieved. Exiting..." + "*".repeat(10));
					flag = false;
					flag2 = false;

					int y = 0;
					for (Bus b : buses) {
						if (b.type == 2) {
							System.out.printf("%.5f \t", b.mp);
							initialPosition[y++] = b.mp;
							
						}
					}
					for (Bus b : buses) {
						if (b.type == 2) {
							System.out.printf("%.5f \t", b.nq);
							initialPosition[y++] = b.nq;
						}
					
					}
					
					System.out.println("\nFull convergence achieved. Exiting..." + "*".repeat(10));
					
				} else if (prev_Mismatches + 0.5e1 <= ModifiedNR.sumMatrix(fx0)) {
					System.out.println("-".repeat(100));

					X1 = prevX1;
					flag = false;
				} else {

					prevX1 = X1;
					prevfx0 = fx0;
					prev_Mismatches = ModifiedNR.sumMatrix(fx0);

					wi = X1.getEntry(X1.getRowDimension() - 2, 0);

					buses.get(0).voltage = new DerivativeStructure(params, order, 1,
							X1.getEntry(X1.getRowDimension() - 1, 0));

					ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
				}

			}

		}

		return initialPosition;

	}
    
    public double[] construct() {

		double[][] xadmittances = this.xadmittances;
		double[][] radmittances = this.radmittances;

		ArrayList<Double> PQLossLoad = null;

		RealMatrix X1 = null;
		RealMatrix fx0 = null;
		boolean flag2 = true;
		double wi = 1.0;
		RealMatrix prevX1 = null;
		RealMatrix prevfx0 = null;
		double prev_Mismatches = Double.MAX_VALUE;
		Random random = new Random();
		boolean flag = true;
		ArrayList<Bus> buses = new ArrayList<Bus>();

		double[][] Jacobian = null;
		double[] initialPosition = null;
		while (flag2) {
			int params = 12;
			int order = 2;

			X1 = null;
			fx0 = null;
			wi = 1.0;
			prevfx0 = null;
			prevX1 = new Array2DRowRealMatrix(params, 1);
			prev_Mismatches = Double.MAX_VALUE;
			flag = true;

			Jacobian = null;
			
			buses = new ArrayList<Bus>();
			buses.add(new Bus(1, 0, params, -0.1 ));
			buses.add(new Bus(1, 2, params, 0.01));
			buses.add(new Bus(1, 4, params, -0.1));
			buses.add(new Bus(2, 6, params, 0, 0, (random.nextDouble() * 2) - 1, (random.nextDouble() * 3), 0.9,12));
			buses.add(new Bus(2, 8, params, 0, 0, (random.nextDouble() * 2) - 1, (random.nextDouble() * 3), 0.9,12));
			buses.add(new Bus(2, 10, params, 0, 0, (random.nextDouble() * 2) - 1, (random.nextDouble() * 3), 0.9,12));

			int nofDroop = 0;
			for (Bus b : buses) {
				if (b.type == 2)
					nofDroop++;
			}

			initialPosition = new double[nofDroop * 2];

			ArrayList<ArrayList<Integer[]>> deltaVoltageOrders = ModifiedNR.createOrders2(buses);
			ArrayList<ArrayList<Integer>> indexes = ModifiedNR.identifyNet(buses);

			innerloop: for (int i = 1; i <= 200 && flag; i++) {
				double w0 = 1.00;
				double v0 = 1.01;

				Complex[][] cAdmittances = Admittance.constructComplexAdmittanceMatrix2(radmittances, xadmittances, wi);

				ModifiedNR.setActiveReactiveGen(buses, wi, w0, v0);

				PQLossLoad = ModifiedNR.calculatePQLossLoad(cAdmittances, buses, wi);

				RealMatrix X0 = ModifiedNR.setUnknownValues(buses, deltaVoltageOrders, wi,
						buses.get(0).voltage.getValue());

				RealMatrix mAdmittances = Admittance.createMadmittance(cAdmittances);

				RealMatrix tAdmittances = Admittance.createTadmittance(cAdmittances);

				ArrayList<ArrayList<DerivativeStructure>> pq = ModifiedNR.createEquations5(buses, mAdmittances,
						tAdmittances);

				double[] mismatches = ModifiedNR.calculateMismatchMatrix2(buses, pq, PQLossLoad, indexes);

				fx0 = new Array2DRowRealMatrix(mismatches);

				Jacobian = ModifiedNR.constructJacobian3(deltaVoltageOrders, pq, buses, wi, cAdmittances, indexes,
						Admittance.createMadmittance(cAdmittances), Admittance.createTadmittance(cAdmittances),
						radmittances, xadmittances);

				RealMatrix JJ = new Array2DRowRealMatrix(Jacobian);

				try {
					try {

						RealMatrix DELTA = MatrixUtils.inverse(JJ).multiply(fx0).scalarMultiply(-1);
						double[][] arJacob = ModifiedNR.artificialJacob(deltaVoltageOrders, buses, wi, indexes,
								radmittances, xadmittances, DELTA, X0);

						RealMatrix aj = new Array2DRowRealMatrix(arJacob);
						DELTA = MatrixUtils.inverse(aj).multiply(fx0).scalarMultiply(-1);
						X1 = X0.add(DELTA);

						if (Double.isNaN(X1.getEntry(X1.getRowDimension() - 2, 0))) {
							System.out.printf("\nAlgorithm diverged, check the line data or the droop coefficients.");
							break innerloop;
						}
					} catch (MatrixDimensionMismatchException | DimensionMismatchException | NullArgumentException e) {

						e.printStackTrace();
					}
				} catch (org.apache.commons.math3.linear.SingularMatrixException e) {

					for (int f = 0; f < Jacobian.length; f++) {
						for (int j = 0; j < Jacobian[1].length; j++) {
							System.out.printf("%.2f \t", Jacobian[f][j]);
						}
						System.out.println();
					}
					for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
						for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
							System.out.printf("%.2f \t", Admittance.createMadmittance(cAdmittances).getEntry(f, j));
						}
						System.out.println();
					}

					for (int f = 0; f < Admittance.createMadmittance(cAdmittances).getRowDimension(); f++) {
						for (int j = 0; j < Admittance.createMadmittance(cAdmittances).getColumnDimension(); j++) {
							System.out.printf("%.2f,%.2f \t", cAdmittances[f][j].getReal(),
									cAdmittances[f][j].getImaginary());
						}
						System.out.println();
					}

					break innerloop;

				}
				if (ModifiedNR.sumMatrix(fx0) <= 1E-4) {

					System.out.println("\nFull convergence achieved. Exiting..." + "*");
					flag = false;
					flag2 = false;

					int y = 0;
					for (Bus b : buses)
						if (b.type == 2)
							initialPosition[y++] = b.mp;

					for (Bus b : buses)
						if (b.type == 2)
							initialPosition[y++] = b.nq;

				} else if (prev_Mismatches + 0.5e2 <= ModifiedNR.sumMatrix(fx0)) {

					System.out.printf("\nMissed the local optimum. Exiting... \t At iteration : %d", i);
					System.out.printf("\nRejected sum of mismatches: %.5f", ModifiedNR.sumMatrix(fx0));
					X1 = prevX1;
					flag = false;
				} else {

					prevX1 = X1;
					prevfx0 = fx0;
					prev_Mismatches = ModifiedNR.sumMatrix(fx0);

					wi = X1.getEntry(X1.getRowDimension() - 2, 0);

					buses.get(0).voltage = new DerivativeStructure(params, order, 1,
							X1.getEntry(X1.getRowDimension() - 1, 0));

					ModifiedNR.updateUnknowns(X1, buses, deltaVoltageOrders, params, order);
				}

			}

		}

		return initialPosition;

	}

    
    /* Sets the conflicts of the honey.
	 *
	 * @param: new number of conflicts
	 */
    public void setConflicts(int mConflicts) {
        this.conflicts = mConflicts;
    }

    /* Gets the selection probability of the honey.
	 *
	 * @return: selection probability of the honey
	 */
    public double getSelectionProbability() {
        return selectionProbability;
    }

    /* sets the selection probability of the honey.
	 *
	 * @param: new selection probability of the honey
	 */
    public void setSelectionProbability(double mSelectionProbability) {
        this.selectionProbability = mSelectionProbability;
    }

    /* Gets the fitness of a honey.
	 *
	 * @return: fitness of honey
	 */
    public double getFitness() {
        return fitness;
    }

    /* Sets the fitness of the honey.
	 *
	 * @param: new fitness
	 */
    public void setFitness(double mFitness) {
        this.fitness = mFitness;
    }

    /* Gets the data on a specified index.
	 *
	 * @param: index of data
	 * @return: position of queen
	 */
    public double getNectar(int index) {
        return nectar[index];
    }

    /* Gets the index on a specified data.
	 *
	 * @param: index of data
	 * @return: position of queen
	 */
    public int getIndex(double value) {
        int k = 0;
        for(; k < MAX_LENGTH; k++) {
            if(nectar[k] == value) {
                break;
            }
        }
        return k;
    }

    /* Sets the data on a specified index.
	 *
	 * @param: index of data
	 * @param: new position of queen
	 */
    public void setNectar(int index, double value) {
        this.nectar[index] = value;
    }

	/* Gets the number of trials of a solution.
	 *
	 * @return: number of trials
	 */
    public int getTrials() {
        return trials;
    }

    /* Sets the number of trials of a solution.
	 *
	 * @param: new number of trials
	 */
    public void setTrials(int trials) {
        this.trials = trials;
    }   
    
    /* Gets the max length.
	 *
	 * @return: max length
	 */
    public int getMaxLength() {
    	return MAX_LENGTH;
    }
}
//end honey