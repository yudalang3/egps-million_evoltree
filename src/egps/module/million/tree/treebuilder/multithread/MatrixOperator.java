package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.multithread;

import java.util.Arrays;
import java.util.Comparator;

import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.TreeBuilder4SarsCov2;
import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit.TreeNodeUnitPair;
import module.remnant.treeoperator.NodeEGPSv1;

public class MatrixOperator {

	static boolean debug;

	private long startTime;

	final private Comparator<TreeNodeUnitPair> pairComparator;

	public MatrixOperator() {
		this(new Comparator<TreeNodeUnitPair>() {
			@Override
			public int compare(TreeNodeUnitPair o1, TreeNodeUnitPair o2) {
				return 0;
			}
		});
	}

	public MatrixOperator(Comparator<TreeNodeUnitPair> pairComparator) {
		this.pairComparator = pairComparator;
		debug = TreeBuilder4SarsCov2.debug;
	}

	public static float[][] doubleArray2float(double[][] dist) {
		int len = dist.length;
		float[][] ret = new float[len][];
		for (int i = 0; i < len; i++) {
			double[] ds = dist[i];
			int len2 = ds.length;
			float[] tt = new float[len2];

			for (int j = 0; j < len2; j++) {
				tt[j] = (float) ds[j];
			}
			ret[i] = tt;
		}

		return ret;
	}

	private class RunningThread extends Thread {
		boolean notStop = true;
		volatile boolean shouldProcess = false;

		int zeroBasedStart = 0;
		int zeroBasedEnd = 10;

		/**
		 * <p>
		 * Title: run
		 * </p>
		 * <p>
		 * Description:
		 * 
		 * A Thread to running background!
		 * 
		 * </p>
		 * 
		 * @see java.lang.Thread#run()
		 *
		 */
		@Override
		public void run() {

			while (notStop) {
				if (shouldProcess) {

					for (int i = zeroBasedStart; i < zeroBasedEnd; i++) {
						for (int j = 0; j < i; j++) {
							if (i > j) {
								//
								float f = distMatrix[i][j];
								if (f<=0) {
									f = 1.1f;
								}
								float ret = (float) (Math.sqrt(f) + Math.pow(f, 3) / 1000 + f*f / f / 4 + f / 30);
//								float ret = f + 2.333f;
								//
								distMatrix[i][j] = ret;
							}

						}
					}

					shouldProcess = false;
				}
			}
		}

		public void operateThread() {
			this.shouldProcess = true;
		}

		public void stopTheThread() {
			notStop = false;
		}

		/**
		 * 
		 * @Title: setOperationRegion @Description: Set region, 0 based!
		 * 
		 * @param: @param start include! @param: @param end exclude! @return:
		 *                void @throws
		 */
		public void setOperationRegion(int start, int end) {
			this.zeroBasedStart = start;
			this.zeroBasedEnd = end;
		}

		public boolean isFinished() {
			return !shouldProcess;
		}

	}

	volatile float[][] distMatrix;
	int lengthOfDistMatrix = 0;
	
	int numOfThread = 3;

	/**
	 * @Title: construct @Description: Main method to build the tree!
	 * 
	 * @param: @param oriTree @param: @param dist @return: void @throws
	 */
	public void construct(NodeEGPSv1 oriTree, float[][] dist) {
		
		RunningThread[] runningThreads = new RunningThread[numOfThread];
		for (int i = 0; i < numOfThread; i++) {
			runningThreads[i] = new RunningThread();
		}

		this.distMatrix = dist;
		lengthOfDistMatrix = dist.length;

		int interval = lengthOfDistMatrix / numOfThread  + 1;
		int start = 0;int end = 0;
		for (int i = 0; i < numOfThread; i++) {
			start = 0 + i * interval;
			end = start + interval;
			if (end > lengthOfDistMatrix) {
				end = lengthOfDistMatrix;
			}
			runningThreads[i].setOperationRegion(start, end);
			runningThreads[i].start();
		}
		
		if (debug) {
			System.out.println("Start!");
		}
		for (int i = 0; i < 30; i++) {

			for (int k = 0; k < numOfThread; k++) {	
				runningThreads[k].operateThread();
			}

			while (true) {
				boolean allFinished = true;
				for (int k = 0; k < numOfThread; k++) {
					if(!runningThreads[k].isFinished()) {
						allFinished = false;
						break;
					}
				}

				// System.out.println("States of thread1 and 2: " + finished1 + "\t" +
				// finished2);
				if (allFinished) {
					break;
				}
			}

			if (debug) {
				System.out.println("Round " + i + " --------------------------");
				System.out.println(Arrays.toString(distMatrix[0]));
			}
			
//			for (float[] fs : distMatrix) {
//				System.out.println(Arrays.toString(fs));
//			}
			

//			
//			try {
//				TimeUnit.SECONDS.sleep(2);
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}

		}

		for (int k = 0; k < numOfThread; k++) {
			runningThreads[k].stopTheThread();
		}
	}
	
	/**  
	 * @Title:  setNumOfThread <BR>  
	 * @Description: please write your description <BR>  
	 * @return: int <BR>  
	 */
	public void setNumOfThread(int numOfThread) {
		this.numOfThread = numOfThread;
	}

	/**
	 * @Title: getDistMatrix <BR>
	 * @Description: please write your description <BR>
	 * @return: float[][] <BR>
	 */
	public float[][] getDistMatrix() {
		return distMatrix;
	}

//    public static void main(String[] args) {
//    	int n = 10;
//		int length = n - 1;
//	    double[][] distance = new double[length][];
//	    
//	    double[] temp0 = {0.0516d};
//	    double[] temp1 = {0.0550d, 0.0031d};
//	    double[] temp2 = {0.0483d, 0.0221d, 0.0253d};
//	    double[] temp3 = {0.0582d, 0.0651d, 0.0685d, 0.0549d};
//	    double[] temp4 = {0.0094d, 0.0416d, 0.0450d, 0.0384d, 0.0549d};
//	    double[] temp5 = {0.0125d, 0.0584d, 0.0619d, 0.0551d, 0.0651d, 0.0157d};
//	    double[] temp6 = {0.0284d, 0.0687d, 0.0722d, 0.0654d, 0.0754d, 0.0317d, 0.0285d};
//	    double[] temp7 = {0.0925d, 0.1221d, 0.1259d, 0.1185d, 0.1370d, 0.0820d, 0.0786d, 0.0927d};
//	    double[] temp8 = {0.1921d, 0.2183d, 0.2228d, 0.2054d, 0.2309d, 0.1798d, 0.1795d, 0.1833d, 0.1860d};
//	    distance[0] = temp0;
//	    distance[1] = temp1;
//	    distance[2] = temp2;
//	    distance[3] = temp3;
//	    distance[4] = temp4;
//	    distance[5] = temp5;
//	    distance[6] = temp6;
//	    distance[7] = temp7;
//	    distance[8] = temp8;
//	    
//	    String[] names = {"seq01", "seq02", "seq03", "seq04", "seq05", "seq06", "seq07", "seq08", "seq09", "seq10"};
//        NJ nj = new NJ();
//        Node root = nj.tree(distance, names);
//        
//        TreeUtility util = new TreeUtility();
//        // root = util.rootAtMidPoint(root);
//        
//        root = util.setRootAt(root, "seq10");
//        
//        System.out.println(egps.module.phylogenetictree.io.TreeCoder.code(root));
//    }

//    public static void main(String[] args) {
//		int n = 11;
//		int length = n - 1;
//		double[][] distance = new double[length][];
//
//		double[] temp1  = { 0.00000d };
//		double[] temp2  = { 0.00299d, 0.00299d };
//		double[] temp3  = { 0.02395d, 0.02395d, 0.02096d };
//		double[] temp4  = { 0.02395d, 0.02395d, 0.02096d, 0.00000d };
//		double[] temp5  = { 0.04491d, 0.04491d, 0.04192d, 0.03892d, 0.03892d };
//		double[] temp6  = { 0.11677d, 0.11677d, 0.11377d, 0.11078d, 0.11078d, 0.07784d };
//		double[] temp7  = { 0.05988d, 0.05988d, 0.05689d, 0.05389d, 0.05389d, 0.01497d, 0.07485d };
//		double[] temp8  = { 0.05988d, 0.05988d, 0.05689d, 0.05389d, 0.05389d, 0.01497d, 0.07485d, 0.00000d };
//		double[] temp9  = { 0.06587d, 0.06587d, 0.06287d, 0.05988d, 0.05988d, 0.03293d, 0.08982d, 0.02994d, 0.02994d };
//		double[] temp10 = { 0.18563d, 0.18563d, 0.18263d, 0.17365d, 0.17365d, 0.15269d, 0.15863d, 0.15269d, 0.15269d, 0.15868d };
//		
//		distance[0 ] = temp1;
//		distance[1 ] = temp2;
//		distance[2 ] = temp3;
//		distance[3 ] = temp4;
//		distance[4 ] = temp5;
//		distance[5 ] = temp6;
//		distance[6 ] = temp7;
//		distance[7 ] = temp8;
//		distance[8 ] = temp9;
//		distance[9 ] = temp10;
//
//		String[] names = { "seq01", "seq02", "seq03", "seq04", "seq05", "seq06", "seq07", "seq08", "seq09", "seq10", "seq11" };
//		NJ nj = new NJ();
//		Node root = nj.tree(distance, names);
//
//		System.out.println(egps.module.phylogenetictree.io.TreeCoder.code(root));
//	}

//  public static void main(String[] args) {
//	  int n = 8;
//		int length = n - 1;
//		double[][] distance = new double[length][];
//
//		double[] temp1 = {1.1639227053989822};
//		double[] temp2 = {1.2450749734986863, 1.2043487537135298};
//		double[] temp3 = {1.14084361917565, 1.250296750156839, 1.3904203140078775};
//		double[] temp4 = {1.0217273959336455, 1.4145031102324044, 1.4650506155883536, 0.7976295040200498};
//		double[] temp5 = {1.1225336130640227, 1.030116700601212, 0.9142262284257272, 1.1793816524372298, 1.39255299900209};
//		double[] temp6 = {1.2462684921299854, 0.6883609452908238, 0.9509558346197838, 1.1782677015451197, 1.2074531793321126, 0.7863096387776374};
//		double[] temp7 = {0.6370331765361604, 1.0175456603655237, 1.676838634880002, 1.0903896177898802, 1.0085478993545594, 1.2952962544491755, 1.1944023318295123};
//		
//		distance[0 ] = temp1;
//		distance[1 ] = temp2;
//		distance[2 ] = temp3;
//		distance[3 ] = temp4;
//		distance[4 ] = temp5;
//		distance[5 ] = temp6;
//		distance[6 ] = temp7;
//
//		String[] names = { "seq01", "seq02", "seq03", "seq04", "seq05", "seq06", "seq07", "seq08" };
//		NJ nj = new NJ();
//		Node root = nj.tree(distance, names);
//		System.out.println(egps.module.phylogenetictree.io.TreeCoder.code(root));
//	}

	public static void main(String[] args) {
		
		int numOfTotalNodes = 1000;
		float[][] dist = new float[numOfTotalNodes][];

		for (int i = 0; i < numOfTotalNodes; i++) {
			float[] tt = new float[i + 1];
			dist[i] = tt;
		}

		MatrixOperator treeAppender = new MatrixOperator();
		treeAppender.construct(null, dist);
	}
}
