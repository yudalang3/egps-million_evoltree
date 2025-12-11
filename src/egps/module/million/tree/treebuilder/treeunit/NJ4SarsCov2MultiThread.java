package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.multithread.RowLocationAssigner;
import evoltree.phylogeny.DefaultPhyNode;
import evoltree.struct.EvolNode;
import evoltree.struct.TreeCoder;
import utils.EGPSUtil;

public class NJ4SarsCov2MultiThread {

	public static boolean debug;
	public static boolean timeTest;

	final private Comparator<TreeNodeUnitPair> pairComparator;

	private volatile float[][] matrixs;
	private volatile List<TreeNodeUnit> listOfNewTreeUnits;
	private volatile double theRestNumOfOTU;
	/** Total number of OTUs + 1! uses[0] is useless; uses[1] is 0.0 */
	private volatile boolean[] uses;
	private volatile int[] nodeIdIndexesOfMatrixTables;
	private volatile double[] sumOfColumns;

	private int numOfThreadIndex = 3;

	private String[] OTUNames;

	private class RunningThread extends Thread {
		volatile boolean notStop = true;
		volatile boolean shouldProcess = false;
		volatile TreeNodeUnitPair minsameQValuePair;
		volatile double tempMinInner;

		int zeroBasedStart = 0;
		int zeroBasedEnd = 10;

		final List<RunningThread> alllRunningThreads;

		public RunningThread(List<RunningThread> alllRunningThreads) {
			this.alllRunningThreads = alllRunningThreads;
		}

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

					long startTime = System.currentTimeMillis();
					minsameQValuePair = null;
					runIterationProcess();
					shouldProcess = false;
					if (!notStop) {
						alllRunningThreads.remove(this);
						// System.err.println("!!!!!!!!!!!!!!!!Thread: " + zeroBasedStart +
						// "Finished!");
						break;
					}
					
					if (timeTest) {
						long endTime = System.currentTimeMillis();
						System.out.println("In Thread " + zeroBasedStart + " Takes: " + (endTime - startTime)
								+ " ms to get min value!");
					}
				}
			}
		}

		/**
		 * @Title: runIterationProcess @Description: 一个线程里面正式做的事情
		 * 
		 * @param: @return: void @throws
		 */
		private void runIterationProcess() {

			/**
			 * 这里就是我们要重点改的地方了，要十分注意各种情况！
			 */
			tempMinInner = Double.MAX_VALUE;
			// mini = minj = 0; // added by fu
			/* Compute Sij values and find the smallest one */
			double firstDivider = theRestNumOfOTU - 2.0;

			List<TreeNodeUnitPair> listOfJointsEachThread = null;
			boolean shouldTeminate = true;
			for (int i = zeroBasedStart; i <= zeroBasedEnd; i++) {
				if (uses[i]) {
					shouldTeminate = false;
					for (int j = 1; j < i; j++) {
						if (uses[j]) {
							/**
							 * 下面这个就是Q value!
							 */
							double SS = firstDivider * matrixs[i][j] - sumOfColumns[j] - sumOfColumns[i];
							
//							double abs = Math.abs(SS - tempMinInner);
//							if (abs > 0 && abs < 0.00001) {
//								System.out.println();
//							}

							if (SS < tempMinInner) {
								TreeNodeUnitPair minsameQValuePair = new TreeNodeUnitPair();
								tempMinInner = SS;
								minsameQValuePair.minI = i;
								minsameQValuePair.minJ = j;
								minsameQValuePair.pairNodeUnitA = listOfNewTreeUnits
										.get(nodeIdIndexesOfMatrixTables[j] - 1);
								minsameQValuePair.pairNodeUnitB = listOfNewTreeUnits
										.get(nodeIdIndexesOfMatrixTables[i] - 1);

								listOfJointsEachThread = new LinkedList<>();
								listOfJointsEachThread.add(minsameQValuePair);
							} else if (SS == tempMinInner) {
								TreeNodeUnitPair other = new TreeNodeUnitPair();
								other.minI = i;
								other.minJ = j;
								other.pairNodeUnitA = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTables[j] - 1);
								other.pairNodeUnitB = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTables[i] - 1);

								listOfJointsEachThread.add(other);

							}
							

						}
					}
				}

			}
			
			if (shouldTeminate) {
				notStop = false;
			}else {
				Iterator<TreeNodeUnitPair> iterator = listOfJointsEachThread.iterator();
				minsameQValuePair = iterator.next();

				while (iterator.hasNext()) {
					TreeNodeUnitPair treeNodeUnitPair = iterator.next();
					if (pairComparator.compare(treeNodeUnitPair, minsameQValuePair) > 0) {
						minsameQValuePair = treeNodeUnitPair;
					}
				}
				
//				int size = listOfJointsEachThread.size();
//				System.out.println(size);
//				if (size < 8) {
//					System.out.println(listOfJointsEachThread);
//				}else {
//					System.out.println(listOfJointsEachThread.subList(0, 3));
//				}
				
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

			if (zeroBasedStart < 2) {
				zeroBasedStart = 2;
			}
		}

		public boolean isFinished() {
			return !shouldProcess;
		}

		/**
		 * @Title: getBestPair @Description: 得到一个线程的最终结果！
		 * 
		 * @param: @return @return: TreeNodeUnitPair @throws
		 */
		public TreeNodeUnitPair getBestPair() {
			return minsameQValuePair;
		}

		public double getTempMinQValue() {
			return tempMinInner;
		}

	}

	public NJ4SarsCov2MultiThread() {
		this(new Comparator<TreeNodeUnitPair>() {
			@Override
			public int compare(TreeNodeUnitPair o1, TreeNodeUnitPair o2) {
				return 0;
			}
		});
	}

	public NJ4SarsCov2MultiThread(Comparator<TreeNodeUnitPair> pairComparator) {
		this.pairComparator = pairComparator;
	}

	private float getValue(float[][] mat, int i, int j) {
		if (i > j) {
			return mat[i][j];
		} else {
			return mat[j][i];
		}

	}

	public void setNumOfThreadIndex(int numOfThreadIndex) {
		this.numOfThreadIndex = numOfThreadIndex;
	}

	/**
	 * Distance matrix should be of the format x1 x2 x3 x4 x5 x6
	 *
	 * where x1 is the distance between sequence 1 and 2
	 * 
	 * @throws IllegalAccessException
	 * @throws InstantiationException
	 */

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

	public static double[][] floatArray2double(float[][] dist) {
		int len = dist.length;
		double[][] ret = new double[len][];
		for (int i = 0; i < len; i++) {
			float[] ds = dist[i];
			int len2 = ds.length;
			double[] tt = new double[len2];

			for (int j = 0; j < len2; j++) {
				tt[j] = ds[j];
			}
			ret[i] = tt;
		}

		return ret;
	}

	// The case tested is from https://en.wikipedia.org/wiki/Neighbor_joining

	public static void main(String[] args) throws ReflectiveOperationException {
		int n = 5;
		int length = n - 1;
		double[][] distance = new double[length][];

		double[] temp0 = { 5.0d };
		double[] temp1 = { 9.0d, 10.0d };
		double[] temp2 = { 9.0d, 10.0d, 8.0d };
		double[] temp3 = { 8.0d, 9.0d, 7.0d, 3.0d };
		distance[0] = temp0;
		distance[1] = temp1;
		distance[2] = temp2;
		distance[3] = temp3;

		String[] names = { "a", "b", "c", "d", "e" };
		NJ4SarsCov2MultiThread nj = new NJ4SarsCov2MultiThread();
		nj.debug = true;
		nj.setNumOfThreadIndex(2);
		TreeNodeUnit tt = new DefaultTreeNodeUnit();
		TreeNodeUnit[] units = { tt, tt, tt, tt, tt };
		EvolNode root = nj.tree(NJ4SarsCov2MultiThread.doubleArray2float(distance), names, units);

		System.out.println();
		String treeString = new TreeCoder<EvolNode>().code(root);
		System.out.println(treeString);

		String realString = "(((a:2.0,b:3.0):3.0,c:4.0):2.0,d:2.0,e:1.0):0.0";
		boolean equals = realString.equals(treeString);

		if (!equals) {
			System.err.println("Error: ture tree is");
			System.err.println(realString);
		}

	}

//	public static void main(String[] args) {
//		int n = 10;
//		int length = n - 1;
//		double[][] distance = new double[length][];
//
//		double[] temp0 = { 0.0516d };
//		double[] temp1 = { 0.0550d, 0.0031d };
//		double[] temp2 = { 0.0483d, 0.0221d, 0.0253d };
//		double[] temp3 = { 0.0582d, 0.0651d, 0.0685d, 0.0549d };
//		double[] temp4 = { 0.0094d, 0.0416d, 0.0450d, 0.0384d, 0.0549d };
//		double[] temp5 = { 0.0125d, 0.0584d, 0.0619d, 0.0551d, 0.0651d, 0.0157d };
//		double[] temp6 = { 0.0284d, 0.0687d, 0.0722d, 0.0654d, 0.0754d, 0.0317d, 0.0285d };
//		double[] temp7 = { 0.0925d, 0.1221d, 0.1259d, 0.1185d, 0.1370d, 0.0820d, 0.0786d, 0.0927d };
//		double[] temp8 = { 0.1921d, 0.2183d, 0.2228d, 0.2054d, 0.2309d, 0.1798d, 0.1795d, 0.1833d, 0.1860d };
//		distance[0] = temp0;
//		distance[1] = temp1;
//		distance[2] = temp2;
//		distance[3] = temp3;
//		distance[4] = temp4;
//		distance[5] = temp5;
//		distance[6] = temp6;
//		distance[7] = temp7;
//		distance[8] = temp8;
//
//		TreeNodeUnit tt = new DefaultTreeNodeUnit();
//		String[] names = { "seq01", "seq02", "seq03", "seq04", "seq05", "seq06", "seq07", "seq08", "seq09", "seq10" };
//		NJ4SarsCov2MultiThread nj = new NJ4SarsCov2MultiThread();
//		NJ4SarsCov2MultiThread.debug = true;
//
//		TreeNodeUnit[] units = { tt, tt, tt, tt, tt, tt, tt, tt, tt, tt };
//		EvolNode root = nj.tree(NJ4SarsCov2.doubleArray2float(distance), names, units);
//
//		String treeString = new TreeCoder().code(root);
//
//		System.out.println(treeString);
//
//		String targetTree = "(seq01:0.0051925815641880035,((((seq02:0.0,seq03:0.0033035704400390387):0.013818749226629734,seq04:0.008331249468028545):0.013251250609755516,seq05:0.03574874997138977):0.015782423317432404,seq06:0.0016675770748406649):0.0028199180960655212,(seq07:0.0059401048347353935,(seq08:0.016810152679681778,(seq09:0.04532499983906746,seq10:0.14067500829696655):0.02818984165787697):0.00383490021340549):0.004344921559095383):0.0";
//		if (!targetTree.equals(treeString)) {
//			System.err.println("Error!  " + "True tree is :");
//			System.err.println(targetTree);
//		}
//	}

	/**
	 * @Title: tree @Description: 重要建树方法的入口
	 * 
	 * @param: @param doubleArray2float @param: @param names @param: @param
	 *                units @param: @return @return: NodeEGPSv1 @throws
	 */
	public EvolNode tree(float[][] dist, String[] OTUnames, TreeNodeUnit[] treeNodeUnits) {
		long startTime = System.currentTimeMillis();

		/**
		 * maximum number of OTU 为什么这里是+1，out不是固定是 +1吗 ？ 那是因为这个程序是原来从一个 1-based
		 * 的计算机程序改编过来的。 因为像 C/C++ Java都是0起始的，但是新的语言一般都是1起始的。 所以这里加了2可以，模拟1起始。
		 */
		final int maxOTU = dist.length + 2;
		/**
		 * parent OTU of each node 这里为什么乘以2，因为二叉树的所有节点数量是叶子的2倍 现在不需要了，parent的值是另一个值的下标。
		 */
		// int[] parent = new int[maxOTU * 2];
		/** branch length of each node */
		float[] branch = new float[maxOTU * 2];
		this.OTUNames = OTUnames;

		List<EvolNode> newTree = new ArrayList<>(16384);
		listOfNewTreeUnits = new ArrayList<>(16384);
		nodeIdIndexesOfMatrixTables = new int[maxOTU];
		uses = new boolean[maxOTU];

		/**
		 * 这个mean是计算一个值的中间变量！
		 */
		double[] mediatedVarMean = new double[maxOTU]; /* ? */

		/**
		 * read and write OTU names and Global initialization 这个到底和 maxOTU 是什么关系！！！
		 * 
		 * 终于知道什么关系了，见第一个注释！
		 * 
		 */
		int nOTU = dist.length + 1;/* No. of OTU (oeprational taxonomic unit) */
		if (debug) {
			System.out.println("This is the debug information from " + getClass());
			System.out.println("Number of OTUs = " + nOTU);
		}

		/**
		 * 这里为什么要从1开始呢？ 因为貌似有这个maxOTU,新加了一列！！！！！！！
		 */
		for (int iOTU = 1; iOTU <= nOTU; iOTU++) {
			// 这个不需要初始化，ydl
			// mat[iOTU][iOTU] = 0.0; /* initialization of diagonal elements */
			uses[iOTU] = true; /* initialization of use vector */
			nodeIdIndexesOfMatrixTables[iOTU] = iOTU; /* initialization of node_id */
			// 这个不需要初始化，ydl
			// mediatedVarMean[iOTU] = 0.0; /* initialization */
		}

		/**
		 * from integer to real 什么玩意，还不如不注释.FOUT就是每次不是会减少一个otu吗？这个存储的就是当前的还剩下的otu的数量
		 * 又删除了下面的废话！！！
		 */
		// for (int i = 1; i <= nOTU * 2 - 3; ++i)
		// parent[i] = 0; /* initialization; OTU 0 is parent for all nodes */

		/* fu's initilization */
		for (int i = 0; i < nOTU; i++) {
			EvolNode node = new DefaultPhyNode();
			node.setName(OTUnames[i]);

			newTree.add(node);
			listOfNewTreeUnits.add(treeNodeUnits[i]);
		}
		matrixs = new float[maxOTU][];
		for (int jOTU = 1; jOTU <= nOTU; jOTU++) {

			float[] ttArray = new float[jOTU + 1];
			for (int iOTU = 1; iOTU < jOTU; iOTU++) {
				ttArray[iOTU] = dist[jOTU - 2][iOTU - 1];
			}
			matrixs[jOTU] = ttArray;
			// matrix[iOTU][jOTU] = dist[jOTU - 2][iOTU - 1];
			// don't need up triangle!
			// matrix[jOTU][iOTU] = dist[jOTU - 2][iOTU - 1];
			/* symmetrization */
		}

//		if (debug) {
//			printDistMatrix();
//		}
		// 分配好多线程的区段：
		RowLocationAssigner rowLocationAssigner = new RowLocationAssigner();
		List<int[]> allocateRegion = rowLocationAssigner.allocateRegionPrecise(nOTU, numOfThreadIndex);

		int numOfThread = allocateRegion.size();

		if (debug || timeTest) {
			System.out.println("Use " + numOfThread + " threads to build the tree!");

			System.out.println("The allocated regions are: ");
			for (int[] is : allocateRegion) {
				System.out.println(Arrays.toString(is));
			}
		}

		CopyOnWriteArrayList<RunningThread> runningThreads = new CopyOnWriteArrayList<>();
		for (int i = 0; i < numOfThread; i++) {
			runningThreads.add(new RunningThread(runningThreads));
		}

		for (int i = 0; i < numOfThread; i++) {
			int[] js = allocateRegion.get(i);
			RunningThread runningThread = runningThreads.get(i);
			runningThread.setOperationRegion(js[0], js[1]);
			runningThread.start();
		}

		/**
		 * specific to NJ method (transformed from PASCAL program) = nc + nOTU
		 * 这个变量像是每一轮新的节点的id！起始应该是新加节点的序号。 例如第一个是nOTU+1,就是新的节点。
		 */
		int theOrderOfNewCreatedNode = 0;
		sumOfColumns = new double[maxOTU];

		// add by yudalang
		for (int i = 1; i <= nOTU; ++i) {
			double temp = 0;

			if (uses[i]) {
				// left
				for (int j = 1; j <= i; j++) {
					if (uses[j]) {
						temp += matrixs[i][j];
					}
				}
				// down
				for (int j = i + 1; j <= nOTU; j++) {
					if (uses[j]) {
						temp += matrixs[j][i];
					}
				}
			}

			sumOfColumns[i] = temp;
		}

		theRestNumOfOTU = nOTU;
		int numOfMainLoop = nOTU - 3;

		for (int nc = 1; nc <= numOfMainLoop; ++nc) {
			if (timeTest) {
				long endTime = System.currentTimeMillis();
				System.out.println("Round " + (nc + 1) + " , total " + numOfMainLoop + "\n-------------------------");
				System.out.println("Takes: " + (endTime - startTime) + " ms to the first of loop;");
				startTime = endTime;
			}
			theOrderOfNewCreatedNode = nc + nOTU;

			// int mini, minj; /* temporary minimum i & j OTU */
//			for (int i = 1; i <= nOTU; ++i) {
//				double temp = 0;
//
//				if (uses[i]) {
//					// left
//					for (int j = 1; j <= i; j++) {
//						if (uses[j]) {
//							temp += matrixs[i][j];
//						}
//					}
//					// down
//					for (int j = i + 1; j <= nOTU; j++) {
//						if (uses[j]) {
//							temp += matrixs[j][i];
//						}
//					}
//				}
//				
//				sumOfColumns[i] = temp;
//			}

			if (timeTest) {
				long endTime = System.currentTimeMillis();
				System.out.print(". Takes: " + (endTime - startTime) + " ms to cal R;");
				startTime = endTime;
			}

			// 多线程计算！！
			for (RunningThread runningThread : runningThreads) {
				runningThread.operateThread();
			}
			while (true) {
				boolean allFinished = true;
				for (RunningThread runningThread : runningThreads) {
					if (!runningThread.isFinished()) {
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

			Iterator<RunningThread> iterator = runningThreads.iterator();
			//第一个一定有
			RunningThread firstRunningThread = iterator.next();
			TreeNodeUnitPair bestPair = firstRunningThread.getBestPair();
			double minQValue = firstRunningThread.getTempMinQValue();
			
			
			while (iterator.hasNext()) {
				RunningThread runningThread = iterator.next();
				TreeNodeUnitPair bestPair2 = runningThread.getBestPair();
				double tempMinQValue = runningThread.getTempMinQValue();

				if (tempMinQValue < minQValue) {
					bestPair = bestPair2;
					minQValue = tempMinQValue;
				} else if (tempMinQValue == minQValue) {
					if (pairComparator.compare(bestPair2, bestPair) > 0) {
						bestPair = bestPair2;
						minQValue = tempMinQValue;
					}
				}
			}

//			if (debug) {
//				System.out.println("Round " + nc + " --------------------------");
//
//				//System.out.println("Sum of colunms are:");
//				//System.out.println(Arrays.toString(sumOfColumns));
//				// System.out.println("Best pair is : " + bestPair);
//				System.out
//						.println("Best pair is : " + OTUnames[bestPair.minI - 1] + "\t" + OTUnames[bestPair.minJ - 1]);
//				System.out.println("Min Q value is :" + minQValue);
//
//			}

			if (timeTest) {
				long endTime = System.currentTimeMillis();

				System.out.println("In main loop, Take " + (endTime - startTime) + " ms to get min pair!");
				startTime = endTime;
			}
			/**
			 * 注意：为了和原来的次序一样，我们这里要交换一下。 它的i是列，j是行
			 */
			int mini = bestPair.minJ;
			int minj = bestPair.minI;
//			int mini = minsameQValuePair.minI;
//			int minj = minsameQValuePair.minJ;

			/**
			 * if (mini == 1 || minj == 1) { int upBound = OTUnames.length - 1;
			 * 
			 * int ti = nodeIdIndexesOfMatrixTables[mini] - 1; int tj =
			 * nodeIdIndexesOfMatrixTables[minj] - 1;
			 * 
			 * ti = ti > upBound ? upBound : ti; tj = tj > upBound ? upBound : tj; String
			 * name1 = OTUnames[ti]; String name2 = OTUnames[tj];
			 * 
			 * if (name1.equals(MainAnalysis.referenceGenomeID) ||
			 * name2.equals(MainAnalysis.referenceGenomeID)) {
			 * System.out.println("####################################");
			 * System.out.print(name1); System.out.println("\t" + name2);
			 * System.out.println(bestPair); }
			 * 
			 * }
			 */

			/**
			 * Compute branch lengths and print the results
			 * 
			 * 
			 */
			double dmin, dio, djo, bri, brj;
			dmin = getValue(matrixs, mini, minj);
			// 这两个是中间变量
			double firstDivider = theRestNumOfOTU - 2.0;
			dio = (sumOfColumns[mini] - dmin) / firstDivider;
			djo = (sumOfColumns[minj] - dmin) / firstDivider;
			bri = (dmin + dio - djo) * 0.5;
			brj = dmin - bri;
			bri = bri - mediatedVarMean[mini];
			brj = brj - mediatedVarMean[minj];

//			if (debug) {
//				String otu_node_i, otu_node_j;
//
//				if (nodeIdIndexesOfMatrixTables[mini] > nOTU)
//					otu_node_i = "node";
//				else
//					otu_node_i = "OTU ";
//				if (nodeIdIndexesOfMatrixTables[minj] > nOTU)
//					otu_node_j = "node";
//				else
//					otu_node_j = "OTU ";
//
//				System.out.println(otu_node_i + nodeIdIndexesOfMatrixTables[mini] + " (" + bri + ") " + otu_node_j + " "
//						+ nodeIdIndexesOfMatrixTables[minj] + " (" + brj + ")");
//
//				System.out.println(" ------------------------------------------");
//			}

			// !! Don't need!
			// parent[nodeIdIndexesOfMatrixTable[mini]] = theOrderOfNewCreatedNode;
			// parent[nodeIdIndexesOfMatrixTable[minj]] = theOrderOfNewCreatedNode;
			branch[nodeIdIndexesOfMatrixTables[mini]] = (float) bri;
			branch[nodeIdIndexesOfMatrixTables[minj]] = (float) brj;

			/* fu's addition */
			EvolNode node = new DefaultPhyNode();

			for (int kk = 0; kk < 2; kk++) {
				int ii = nodeIdIndexesOfMatrixTables[mini] - 1;
				if (kk == 1)
					ii = nodeIdIndexesOfMatrixTables[minj] - 1;

				EvolNode child = newTree.get(ii);
				node.addChild(child);
			}
			newTree.add(node);

			TreeNodeUnit treeNodeUnitA = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTables[mini] - 1);
			TreeNodeUnit treeNodeUnitB = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTables[minj] - 1);
			listOfNewTreeUnits.add(treeNodeUnitA.produceCluster(treeNodeUnitB));
			// node.setName("OTU " + newTree.size());

			if (timeTest) {
				long endTime = System.currentTimeMillis();

				System.out.print(" Take " + (endTime - startTime) + " ms to cal branch and produce Cluster!");
				startTime = endTime;
			}

			/* re-initialization */
			theRestNumOfOTU = theRestNumOfOTU - 1.0;
			mediatedVarMean[mini] = dmin * 0.5;
			uses[minj] = false;
			// System.out.println("To remove id : "+minj+"
			// "+nodeIdIndexesOfMatrixTable[minj]);
			nodeIdIndexesOfMatrixTables[mini] = theOrderOfNewCreatedNode;

			/* compute new distances */
			double sumOfColumnMinI = 0;
			for (int j = 1; j <= nOTU; ++j) {
				if (uses[j]) {
					float da = getValue(matrixs, mini, j) + getValue(matrixs, minj, j);
					da = 0.5f * da;
					if (mini < j) {
						matrixs[j][mini] = da;
					} else if (mini > j) {
						matrixs[mini][j] = da;
					}
					/**
					 * at the same time update the matrix
					 */
					if (j != mini) {
						sumOfColumns[j] -= da;
						sumOfColumnMinI += da;
					}
				}
			}
			sumOfColumns[mini] = sumOfColumnMinI;
			/**
			 * 被选中的原来的值赋为空 这个不用动啊！！
			 */
//			for (int j = 1; j <= nOTU; ++j) {
//				//matrix[minj][j] = 0f;
//				//matrix[j][minj] = 0f;
//				
//				if (minj < j) {
//					matrix[j][minj] = 0f;
//				}else {
//					matrix[minj][j] = 0f;
//				}
//			}

			if (debug) {
				// printDistMatrix();
				System.out.println(
						"This is round " + (nc + 1) + " , in total main loop, total iteration is " + numOfMainLoop);
			}

			if (timeTest) {
				long endTime = System.currentTimeMillis();

				System.out.println(" Take " + (endTime - startTime) + "ms to last code in main loop!");
				startTime = endTime;
			}

		}
		/* end of the main for loop */
		for (RunningThread runningThread : runningThreads) {
			runningThread.stopTheThread();
		}

		/**
		 * The last cycle
		 */
		int[] L = new int[4]; /* OTU ID of the last 3-OTU tree */
		double[] LB = new double[4]; /* branch lengths of the last 3-OTU tree */

		int jun = 1;
		for (int i = 1; i <= nOTU; ++i) {
			if (uses[i]) {
				L[jun] = i;
				++jun;
			}
		}

		float v12 = getValue(matrixs, L[1], L[2]);
		float v23 = getValue(matrixs, L[2], L[3]);
		float v13 = getValue(matrixs, L[1], L[3]);
		LB[1] = (v12 + v13 - v23) * 0.5;
		LB[2] = v12 - LB[1];
		LB[3] = v13 - LB[1];

		EvolNode last = new DefaultPhyNode();
		// last.setName("OTU " + (newTree.size() + 1));
		for (int i = 1; i <= 3; ++i) {
			LB[i] = LB[i] - mediatedVarMean[L[i]];

//			if (debug) {
//				String otu_node_i;
//				if (nodeIdIndexesOfMatrixTables[L[i]] > nOTU)
//					otu_node_i = "node ";
//				else
//					otu_node_i = "OTU ";
//				System.out.println(otu_node_i + nodeIdIndexesOfMatrixTables[L[i]] + " (" + LB[i] + " )");
//			}

			branch[nodeIdIndexesOfMatrixTables[L[i]]] = (float) LB[i];

			/*
			 * this part is added by Fu Note that the last node is a trifurcation
			 */

			EvolNode child = newTree.get(nodeIdIndexesOfMatrixTables[L[i]] - 1);
			last.addChild(child);
		}
		
		
		newTree.add(last);

		/* end of computation */

		/* print branch lengths */
		// if(debug)System.out.println("\n== Node length list (Node 0 is the last node)
		// ==");
		for (int i = 1; i <= nOTU * 2 - 3; ++i) {
			// if (debug) System.out.println("branch "+i+" ->"+ parent[i]+" ("+branch[i]+"
			// )");
			EvolNode node = newTree.get(i - 1);

			branch[i] = (branch[i] < 0 ? 0 : branch[i]); // Added by Haipeng

			node.setLength(branch[i]);
		}
		
		//added by yudalang :trifurcation should be a bifurcation
		EvolNode firstChild = last.removeChild(0);
		double length = firstChild.getLength();
		double halfLen = 0.5 * length;
		firstChild.setLength(halfLen);
		last.setLength(halfLen);
		
		EvolNode root = new DefaultPhyNode();
		root.addChild(firstChild);
		root.addChild(last);
		
		if (debug) {
			System.out.println("Total process finished!");
			startTime = EGPSUtil.printTimeSinceStartAndPrintUsedMemory(startTime);
		}

		return root;
	}

	private void printDistMatrix() {
		int length = matrixs.length;
		if (length < 30) {
			System.out.println("The distance matrix:");

			for (int i = 1; i < length; i++) {
				if (uses[i]) {
					System.out.print(OTUNames[i - 1] + "  ");
					float[] ds = matrixs[i];
					for (int j = 1; j < ds.length; j++) {
						System.out.print(Float.toString(ds[j]).concat(" , "));
					}
					System.out.println();
				}

			}
		}
		System.out.println("------------------------------");
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
}
