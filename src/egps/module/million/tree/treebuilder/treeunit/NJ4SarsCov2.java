package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.TreeBuilder4SarsCov2;

import evoltree.phylogeny.DefaultPhyNode;
import evoltree.struct.TreeCoder;
import utils.EGPSUtil;

public class NJ4SarsCov2 {

	static boolean debug;

	private long startTime;

	final private Comparator<TreeNodeUnitPair> pairComparator;

	public NJ4SarsCov2() {
		this(new Comparator<TreeNodeUnitPair>() {
			@Override
			public int compare(TreeNodeUnitPair o1, TreeNodeUnitPair o2) {
				return 0;
			}
		});
	}

	public NJ4SarsCov2(Comparator<TreeNodeUnitPair> pairComparator) {
		this.pairComparator = pairComparator;
		debug = TreeBuilder4SarsCov2.debug;
	}

	/**
	 * Distance matrix should be of the format x1 x2 x3 x4 x5 x6
	 *
	 * where x1 is the distance between sequence 1 and 2
	 * 
	 * @throws IllegalAccessException
	 * @throws InstantiationException
	 */
	public DefaultPhyNode tree(float[][] dist, String[] OTUnames, TreeNodeUnit[] treeNodeUnits) {
		startTime = System.currentTimeMillis();

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

		List<DefaultPhyNode> newTree = new ArrayList<>(16384);
		List<TreeNodeUnit> listOfNewTreeUnits = new ArrayList<>(16384);

		/**
		 * node ID starting from nOTU+1
		 * 
		 * 完全不知道这个注释是什么东西： 下面不是有赋值说 node_id[i] = i 吗？为什么这里有什么from nOTU + 1
		 * 
		 * 第二次的感受，这个变量的作用是得到距离矩阵对应下标的叶子节点或者新的cluster的下标，
		 * 
		 * 注意它的值是下标，不是id，是一个数组的下标。
		 * 
		 * 例如一开始全部的OTU都还没有聚集起来，那么第一个叶子就是1
		 * 
		 * 这里的node包括leaf和internal node
		 */
		int[] nodeIdIndexesOfMatrixTable = new int[maxOTU];

		float[][] matrix = new float[maxOTU][maxOTU]; /* distance matrix */

		/**
		 * OTU usage vector 暂时还不知道这是个什么东西，还好下面又有一条注释： 为什么不用boolean，通俗易懂。
		 * 
		 * 0 = do not use OTU; 1 = use OTU
		 * 
		 * 你会看到有些时候会有起始下标为2的操作 这个值是会变的！！！
		 * 
		 * yudalang： 改成了boolen数组！
		 */
		boolean[] use = new boolean[maxOTU];

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
			use[iOTU] = true; /* initialization of use vector */
			nodeIdIndexesOfMatrixTable[iOTU] = iOTU; /* initialization of node_id */
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
            DefaultPhyNode node = new DefaultPhyNode();
			if (OTUnames == null) {
				node.setName("OTU" + (i + 1));
			} else {
				node.setName(OTUnames[i]);
			}

			newTree.add(node);

			listOfNewTreeUnits.add(treeNodeUnits[i]);
		}
		// 这里有对角线的赋值！
		/* read distance matrix */
		for (int jOTU = 2; jOTU <= nOTU; jOTU++) {
			for (int iOTU = 1; iOTU <= jOTU - 1; iOTU++) {
				matrix[iOTU][jOTU] = dist[jOTU - 2][iOTU - 1];
				matrix[jOTU][iOTU] = dist[jOTU - 2][iOTU - 1];
				/* symmetrization */
			}
		}

		if (debug) {
			System.out.println("The distance matrix: Note, first line and column are empty!!");
			if (maxOTU < 30) {
				for (float[] ds : matrix) {
					System.out.println(Arrays.toString(ds));
				}
			}
			System.out.println("Finished preparing ... " + getClass());
			startTime = System.currentTimeMillis();;
			System.out.println(EGPSUtil.getAlreadyUsedJVMMemory() + " MB memory space has already use!");
		}

		/**
		 * specific to NJ method (transformed from PASCAL program) = nc + nOTU
		 * 这个变量像是每一轮新的节点的id！起始应该是新加节点的序号。 例如第一个是nOTU+1,就是新的节点。
		 */
		int theOrderOfNewCreatedNode = 0;
		/**
		 * sum of Dij for all j 简单的说就是列和！
		 */
		double[] R = new double[maxOTU];

		/* main loop */
		double theRestNumOfOTU = nOTU;
		int numOfMainLoop = nOTU - 3;

		double tempmin; /* temporary minimum distance */
		int mini, minj; /* temporary minimum i & j OTU */

		for (int nc = 1; nc <= numOfMainLoop; ++nc) {
			theOrderOfNewCreatedNode = nc + nOTU;
//			if (debug) {
//				System.out.print("Node " + nnode + " :");
//			}

			/**
			 * Local initialization SUMD不需要，直接删掉
			 */
			// SUMD = 0.0;
			for (int j = 2; j <= nOTU; ++j)
				for (int i = 1; i <= j - 1; ++i) {
					matrix[j][i] = matrix[i][j];
					// SUMD = SUMD + mat[i][j];
				}
			/**
			 * 下面这个我替换成了别的最大值，这个值其实我觉得还是有风险的！
			 */
			// tempmin = 99999999.0;
			tempmin = Integer.MAX_VALUE;
			for (int i = 1; i <= nOTU; ++i) {
				R[i] = 0.0;
				for (int j = 1; j <= nOTU; ++j)
					R[i] = R[i] + matrix[i][j];
			}

			//System.out.println("R array are: " + Arrays.toString(R));

			/**
			 * 这里就是我们要重点改的地方了，要十分注意各种情况！
			 */
			TreeNodeUnitPair minsameQValuePair = new TreeNodeUnitPair();
			int numOfSamePairs = 0;
			mini = minj = 0; // added by fu
			/* Compute Sij values and find the smallest one */
			double firstDivider = theRestNumOfOTU - 2.0;
			for (int jOTU = 2; jOTU <= nOTU; jOTU++) {
				if (use[jOTU]) {
					for (int iOTU = 1; iOTU <= jOTU - 1; iOTU++) {
						if (use[iOTU]) {
							/**
							 * 下面这个就是Q value!
							 */
							double SS = firstDivider * matrix[iOTU][jOTU] - R[iOTU] - R[jOTU];
							if (SS < tempmin) {
								tempmin = SS;
								mini = iOTU;
								minj = jOTU;
								minsameQValuePair.minI = iOTU;
								minsameQValuePair.minJ = jOTU;
								minsameQValuePair.pairNodeUnitA = listOfNewTreeUnits
										.get(nodeIdIndexesOfMatrixTable[iOTU] - 1);
								minsameQValuePair.pairNodeUnitB = listOfNewTreeUnits
										.get(nodeIdIndexesOfMatrixTable[jOTU] - 1);

								numOfSamePairs = 0;
							} else if (SS == tempmin) {
								numOfSamePairs++;

								TreeNodeUnitPair other = new TreeNodeUnitPair();
								other.minI = iOTU;
								other.minJ = jOTU;
								other.pairNodeUnitA = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTable[iOTU] - 1);
								other.pairNodeUnitB = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTable[jOTU] - 1);
								if (pairComparator.compare(other, minsameQValuePair) > 0) {
									minsameQValuePair = other;
								}
							}

							// 不需要用这个Sij 直接用s就行了
							// Sij = (SS * 0.5 + SUMD) / (FOTU - 2.0);
							// if (Sij < tempmin) {
							// tempmin = Sij;
							// mini = iOTU;
							// minj = jOTU;
							// }
						}
					}
				}
			}

			if (debug) {
				System.out.println("Same q value, size is : " + numOfSamePairs);
			}
			mini = minsameQValuePair.minI;
			minj = minsameQValuePair.minJ;

			/**
			 * Compute branch lengths and print the results
			 * 
			 * 
			 */
			double dmin, dio, djo, bri, brj, da;
			dmin = matrix[mini][minj];
			// 这两个是中间变量
			dio = (R[mini] - dmin) / firstDivider;
			djo = (R[minj] - dmin) / firstDivider;
			bri = (dmin + dio - djo) * 0.5;
			brj = dmin - bri;
			bri = bri - mediatedVarMean[mini];
			brj = brj - mediatedVarMean[minj];

			if (debug) {
				String otu_node_i, otu_node_j;

				if (nodeIdIndexesOfMatrixTable[mini] > nOTU)
					otu_node_i = "node";
				else
					otu_node_i = "OTU ";
				if (nodeIdIndexesOfMatrixTable[minj] > nOTU)
					otu_node_j = "node";
				else
					otu_node_j = "OTU ";

				System.out.println(otu_node_i + nodeIdIndexesOfMatrixTable[mini] + " (" + bri + ") " + otu_node_j + " "
						+ nodeIdIndexesOfMatrixTable[minj] + " (" + brj + ")");
			}

			// !! Don't need!
			// parent[nodeIdIndexesOfMatrixTable[mini]] = theOrderOfNewCreatedNode;
			// parent[nodeIdIndexesOfMatrixTable[minj]] = theOrderOfNewCreatedNode;
			branch[nodeIdIndexesOfMatrixTable[mini]] = (float) bri;
			branch[nodeIdIndexesOfMatrixTable[minj]] = (float) brj;

			/* fu's addition */
			DefaultPhyNode node = new DefaultPhyNode();

			for (int kk = 0; kk < 2; kk++) {
				int ii = nodeIdIndexesOfMatrixTable[mini] - 1;
				if (kk == 1)
					ii = nodeIdIndexesOfMatrixTable[minj] - 1;

				DefaultPhyNode child = (DefaultPhyNode) newTree.get(ii);
				node.addChild(child);
			}
			newTree.add(node);
			TreeNodeUnit treeNodeUnitA = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTable[mini] - 1);
			TreeNodeUnit treeNodeUnitB = listOfNewTreeUnits.get(nodeIdIndexesOfMatrixTable[minj] - 1);
			listOfNewTreeUnits.add(treeNodeUnitA.produceCluster(treeNodeUnitB));
			node.setName("OTU " + newTree.size());

			/* re-initialization */
			theRestNumOfOTU = theRestNumOfOTU - 1.0;
			mediatedVarMean[mini] = dmin * 0.5;
			use[minj] = false;
			//System.out.println("To remove id : " + minj + "   " + nodeIdIndexesOfMatrixTable[minj]);
			nodeIdIndexesOfMatrixTable[mini] = theOrderOfNewCreatedNode;

			/* compute new distances */
			for (int j = 1; j <= nOTU; ++j) {
				if (use[j]) {
					da = (matrix[mini][j] + matrix[minj][j]) * 0.5;
					if (mini < j)
						matrix[mini][j] = (float) da;
					else if (mini > j)
						matrix[j][mini] = (float) da;

					//System.out.println("da\t" + da);
				}
			}
			/**
			 * 被选中的原来的值赋为空
			 */
			for (int j = 1; j <= nOTU; ++j) {
				matrix[minj][j] = 0f;
				matrix[j][minj] = 0f;
			}

			if (debug) {
				int length = matrix.length;
				if (length < 30) {

					System.out.println("New matrix is:");

					for (float[] ds : matrix) {
						System.out.println(Arrays.toString(ds));
					}
				}

				System.out.println("This is round " + nc + " , in total main loop, total iteration is " + numOfMainLoop);
				startTime = System.currentTimeMillis();
			}

		}
		/* end of the main for loop */

		/**
		 * The last cycle
		 */
		int[] L = new int[4]; /* OTU ID of the last 3-OTU tree */
		double[] LB = new double[4]; /* branch lengths of the last 3-OTU tree */

		int jun = 1;
		for (int i = 1; i <= nOTU; ++i) {
			if (use[i]) {
				L[jun] = i;
				++jun;
			}
		}
		LB[1] = (matrix[L[1]][L[2]] + matrix[L[1]][L[3]] - matrix[L[2]][L[3]]) * 0.5;
		LB[2] = matrix[L[1]][L[2]] - LB[1];
		LB[3] = matrix[L[1]][L[3]] - LB[1];

		DefaultPhyNode last = new DefaultPhyNode();
		last.setName("OTU " + (newTree.size() + 1));
		for (int i = 1; i <= 3; ++i) {
			LB[i] = LB[i] - mediatedVarMean[L[i]];

			if (debug) {
				String otu_node_i;
				if (nodeIdIndexesOfMatrixTable[L[i]] > nOTU)
					otu_node_i = "node ";
				else
					otu_node_i = "OTU ";
				System.out.println(otu_node_i + nodeIdIndexesOfMatrixTable[L[i]] + " (" + LB[i] + " )");
			}

			branch[nodeIdIndexesOfMatrixTable[L[i]]] = (float) LB[i];

			/*
			 * this part is added by Fu Note that the last node is a trifurcation
			 */

			DefaultPhyNode child = (DefaultPhyNode) newTree.get(nodeIdIndexesOfMatrixTable[L[i]] - 1);
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
			DefaultPhyNode node = (DefaultPhyNode) newTree.get(i - 1);

			branch[i] = (branch[i] < 0 ? 0 : branch[i]); // Added by Haipeng

            node.setLength(branch[i]);; // fu's addition
		}

		if (debug) {
			System.out.println("Total process finished!");
			startTime = System.currentTimeMillis();
			System.out.println(EGPSUtil.getAlreadyUsedJVMMemory() + " MB memory space has already use!");
		}

		return last;
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
		NJ4SarsCov2 nj = new NJ4SarsCov2();
		nj.debug = true;
		TreeNodeUnit tt = new DefaultTreeNodeUnit();
		TreeNodeUnit[] units = { tt, tt, tt, tt, tt, tt };
		DefaultPhyNode root = nj.tree(NJ4SarsCov2.doubleArray2float(distance), names, units);

		System.out.println();
        String treeString = new TreeCoder<DefaultPhyNode>().code(root);
		System.out.println(treeString);

		String realString = "(((a:2.0,b:3.0):3.0,c:4.0):2.0,d:2.0,e:1.0):0.0;";
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
//		NJ4SarsCov2 nj = new NJ4SarsCov2();
//
//		TreeNodeUnit[] units = { tt, tt, tt, tt, tt, tt, tt, tt, tt, tt };
//		DefaultPhyNode root = nj.tree(NJ4SarsCov2.doubleArray2float(distance), names, units);
//
//		String treeString = egps.module.phylogenetictree.io.TreeCoder.code(root);
//
//		System.out.println(treeString);
//	}

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
