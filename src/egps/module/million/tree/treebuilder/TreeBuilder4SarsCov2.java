package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import evoltree.phylogeny.DefaultPhyNode;
import evoltree.struct.EvolNode;
import evoltree.struct.TreeCoder;
import evoltree.struct.util.EvolTreeOperator;
import module.evolview.gfamily.work.CGBNodeUtil;
import module.remnant.treeoperator.NodeEGPSv1;
import module.remnant.treeoperator.reconAlgo.NJ;
import org.apache.commons.io.FileUtils;

import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit.NJ4SarsCov2MultiThread;
import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit.SarsCov2TreeUnitPairComparatorSimple;
import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit.TreeNodeUnit;
import utils.EGPSUtil;

/**
 * 
 * 计算最简单的进化距离，然后用NJ法构一棵NJ树。 进化距离的计算，fasta文件可以分批。支持外群。没有外群也能建树，但是不会reroot。
 * NJ法经过优化，支持多线程。
 * 
 * @title TreeBuilder4SarsCov2
 * @createdDate 2020-10-21 16:47
 * @lastModifiedDate 2020-10-21 16:47
 * @author "yudalang"
 * @since 1.7
 *
 */
public class TreeBuilder4SarsCov2 {

	File fastaDir;
	File totalInforFile;
	File referenceGenomeFile;
	File ancesterStateFile;
	File fastaFileContainsOutgroup;

	String refGenomeID;

	/** 默认的进化树输出路径 */
	File outputNwkFile = new File("out.nwk");
	private File newTotalInforFile = new File("newTotalInfor.txt");

	public static boolean debug = false;

	int realStartPosition = 10;
	int realEndPostion = 20;

	long startTime;
	private int numOfThreads = 3;

	public void setFastaDir(File fastaDir) {
		this.fastaDir = fastaDir;
	}

	public void setNewTotalInforFile(File newTotalInforFile) {
		this.newTotalInforFile = newTotalInforFile;
	}

	public void setTotalInforFile(File totalInforFile) {
		this.totalInforFile = totalInforFile;
	}

	/**
	 * @param fastaFileContainsOutgroup the {@link #fastaFileContainsOutgroup} to
	 *                                  set
	 */
	public void setFastaFileContainsOutgroup(File fastaFileContainsOutgroup) {
		this.fastaFileContainsOutgroup = fastaFileContainsOutgroup;
	}

	public void setReferenceGenomeFile(File referenceGenomeFile) {
		this.referenceGenomeFile = referenceGenomeFile;
	}

	public void setAncesterStateFile(File ancesterStateFile) {
		this.ancesterStateFile = ancesterStateFile;
	}

	public void setRefGenomeID(String refGenomeID) {
		this.refGenomeID = refGenomeID;
	}

	public void setRealCalculatePostion(int start, int end) {
		realStartPosition = start;
		realEndPostion = end;
	}

	public void setNumOfThreads(int numOfThreadIndex) {
		this.numOfThreads = numOfThreadIndex;
	}

	public void setOutputNwkFile(File outputNwkFile) {
		this.outputNwkFile = outputNwkFile;
	}

	public void runAnalysis() throws Exception {
		boolean directory = fastaDir.isDirectory();
		if (!directory) {
			throw new IllegalArgumentException("Should input a fasta dir!");
		}
		File[] listFiles = fastaDir.listFiles();
		// 计算进化距离
		AbstractEvolDistWholeRunnerSingleThread wholeRunner = new EvolDistWholeFromScratch();
		wholeRunner.setInputFastaFiles(Arrays.asList(listFiles));
		wholeRunner.setRefGenomeFile(referenceGenomeFile);
		wholeRunner.setRefGenomeID(refGenomeID);
		wholeRunner.setAncesterStateFile(ancesterStateFile);
		wholeRunner.setRealCalculatePostion(realStartPosition, realEndPostion);
		wholeRunner.setFastaFileContainOutgroupFile(fastaFileContainsOutgroup);

		System.out.println("Start to calculate evolutionary distance ...");
		startTime = System.currentTimeMillis();

		float[][] calculateDistanceMatrix = wholeRunner.calculateDistanceMatrix();

		System.out.println("Finished calculating matrix array...");

		startTime = EGPSUtil.printTimeSinceStartAndPrintUsedMemory(startTime);

		String[] otuNames = wholeRunner.getOUTNames();
		Map<String, Integer> seqName2seqNameIndexMap = wholeRunner.getSeqName2seqNameIndexMap();

		debug2seeTheResults(calculateDistanceMatrix, otuNames, seqName2seqNameIndexMap);

		// Next is build tree
		buildTreeDirectly(calculateDistanceMatrix, otuNames);

		EGPSUtil.printTimeSinceStartAndPrintUsedMemory(startTime);
		System.out.println("Last procedure finished, total process finished!");
	}

	private void debug2seeTheResults(float[][] calculateDistanceMatrix, String[] otuNames,
			Map<String, Integer> seqName2seqNameIndexMap) throws IOException {
		if (debug) {
			/**
			 * 可以在这里查看你想要查看的两条序列之间的距离。
			 */

			// 输出外群之间的距离
//			Integer o1 = seqName2seqNameIndexMap.get("Outgroup_RaTG13");
//			Integer o2 = seqName2seqNameIndexMap.get("Outgroup_GX-P1E");
//
//			float tt = 0;
//			if (o2 > o1) {
//				tt = calculateDistanceMatrix[o2][o1];
//			} else {
//				tt = calculateDistanceMatrix[o1][o2];
//			}
//			System.out.println("Distance of two outGroup: " + tt);

			List<String> outStrings = new ArrayList<String>(calculateDistanceMatrix.length);

			// output all distances
//			StringBuilder sBuilder = new StringBuilder(100000);
//			sBuilder.append(",");
//			for (String string : outNames) {
//				sBuilder.append(string).append(",");
//			}
//			outStrings.add(sBuilder.toString());
//
//			int index = 0;
//			for (float[] fs : calculateDistanceMatrix) {
//				sBuilder.setLength(0);
//				sBuilder.append(outNames[index]).append(",");
//				for (float ft : fs) {
//					sBuilder.append(Float.toString(ft)).append(",");
//
//				}
//				outStrings.add(sBuilder.toString());
//				index++;
//			}

			// output target seq distances
			Set<Integer> targetSeqIndexes = new HashSet<>();
			Integer integer1 = seqName2seqNameIndexMap.get("NC_045512");
			Integer integer2 = seqName2seqNameIndexMap.get("EPI_ISL_417390");
			Integer integer3 = seqName2seqNameIndexMap.get("EPI_ISL_416334");
			Integer integer4 = seqName2seqNameIndexMap.get("EPI_ISL_420346");
			targetSeqIndexes.add(integer1);
//			targetSeqIndexes.add(integer2);
//			targetSeqIndexes.add(integer3);
//			targetSeqIndexes.add(integer4);

//			int length3 = calculateDistanceMatrix.length;
//			for (int i = 0; i < length3; i++) {
//				for (int j = 0; j < i; j++) {
//					if (targetSeqIndexes.contains(i) || targetSeqIndexes.contains(j)) {
//						float f = calculateDistanceMatrix[i][j];
//						outStrings.add(outNames[i] + ","+outNames[j] +","+String.valueOf(f));
//					}
//				}
//			}

			// output distance of refGenome to others
//			StringBuilder sBuilder = new StringBuilder(100000);
//			for (String string : outNames) {
//				sBuilder.setLength(0);
//
//				sBuilder.append(string).append(",");
//
//				Integer integer = seqName2seqNameIndexMap.get(string);
//
//				float dd = calculateDistanceMatrix[integer][0];
//				sBuilder.append(String.valueOf(dd)).append(",");
//				sBuilder.append(readVirusStrainsInfo4Builder.getVirusStrain(string).getDateString());
//
//				outStrings.add(sBuilder.toString());
//			}

//			Map<Float, Integer> map = new HashMap<>();
//			int size = calculateDistanceMatrix.length;
//			for (int i = 0; i < size; i++) {
//				for (int j = 0; j < i; j++) {
//					float ft = calculateDistanceMatrix[i][j];
//					Integer count = map.get(ft);
//					map.put(ft, (count == null) ? 1 : count + 1);
//				}
//			}
//
//			for (Entry<Float, Integer> entry : map.entrySet()) {
//				String ss = entry.getKey() + "," + entry.getValue();
//				outStrings.add(ss);
//			}

			FileUtils.writeLines(new File(outputNwkFile.getAbsolutePath() + ".build.dist.csv"), outStrings);
			// throw new IOException();

		}

	}

	private EvolNode buildTreeDirectly(float[][] calculateDistanceMatrix, String[] outNames) throws Exception {

		boolean hasPloymorsim = false;
		for (float[] string : calculateDistanceMatrix) {
			for (float string2 : string) {
				if (string2 != 0) {
					hasPloymorsim = true;
					break;
				}
			}
			if (hasPloymorsim) {
				break;
			}
		}

		if (!hasPloymorsim) {
			System.err.println("No polymorphism in your MSA, please validate your data!");
			System.exit(0);
		}

		// Should let EvolDistWholeRunner become null!
		int length = calculateDistanceMatrix.length;
		int totalOTULength = outNames.length;

		if (totalOTULength != (length)) {
			System.err.println("Sorry, this this a bug! Please tell yudalang! " + getClass());
			throw new InternalError();
		}
		// 在分析之前先读入基本的毒株信息

		calculateDistanceMatrix = Arrays.copyOfRange(calculateDistanceMatrix, 1, totalOTULength);

		System.out.println(EGPSUtil.getAlreadyUsedJVMMemory() + " MB memory space has already use!");

		/**
		 * Start to build the tree
		 */
		// 单线程
//		 NJMethod_UseFloat_YDL nj = new NJMethod_UseFloat_YDL();
		// 多线程
		NJ4SarsCov2MultiThread nj = new NJ4SarsCov2MultiThread(new SarsCov2TreeUnitPairComparatorSimple());
		nj.setNumOfThreadIndex(numOfThreads);
		TreeNodeUnit[] units = new TreeNodeUnit[totalOTULength];
		List<String> outputStringsOfTotalInfor = new ArrayList<String>(totalOTULength);

//		LocalDate earliestCollectionDate = readVirusStrainsInfo4Builder.getEarliestCollectionDate();
//
//		for (int i = 0; i < totalOTULength; i++) {
//			String accessionNumber = outNames[i];
//
//			if (OutgroupOperator.ifMeetOutgroupPrefix(accessionNumber)) {
//				units[i] = new SarsCov2TreeUnitSimple(5000, "NoCountry", new HashMap<String, Double>());
//			} else {
//				VirusStrain virusStrain = readVirusStrainsInfo4Builder.getVirusStrain(accessionNumber);
//				if (virusStrain == null) {
//					System.err.println(accessionNumber + " don not have virus strain inforamtion!");
//					throw new IOException();
//				}
//
//				outputStringsOfTotalInfor.add(virusStrain.oriLine);
//
//				LocalDate date2 = virusStrain.getCollectionDate();
//
//				int days = (int) (ChronoUnit.DAYS.between(earliestCollectionDate, date2));
//
//				Map<String, Double> country2weightMap = new HashMap<>();
//				country2weightMap.put(virusStrain.country, 1.0);
//				units[i] = new SarsCov2TreeUnitSimple(days, virusStrain.country, country2weightMap);
//			}
//		}
//
//		FileUtils.writeLines(newTotalInforFile, outputStringsOfTotalInfor);
//
//		// clear garbage again
//		readVirusStrainsInfo4Builder = null;

		/**
		 * 下面可以是建树
		 */
		EvolNode root = nj.tree(calculateDistanceMatrix, outNames, units);
		// EvolNode root = reconstructTree(calculateDistanceMatrix,outNames);

		System.out.println("Get root node, used ");
		startTime = EGPSUtil.printTimeSinceStartAndPrintUsedMemory(startTime);

		System.out.println();
		TreeCoder<EvolNode> treeCoder = new TreeCoder<EvolNode>();
		String unRerootedTreeString = treeCoder.code(root);

		/**
		 * 重新定根在输出树
		 */
		List<EvolNode> outgroups = new ArrayList<>();
		final String outGroupName = "Outgroup";
		List<EvolNode> leaves = CGBNodeUtil.getLeavesByRecursive(root);
		for (EvolNode evolNode : leaves) {
			if (evolNode.getName().startsWith(outGroupName)) {
				outgroups.add(evolNode);
			}
		}

		EvolNode reNode = null;
		if (outgroups.size() == 0) {
			System.err.println("Out group root error! Cannot find the outgroup node.");
		} else if (outgroups.size() == 1) {
			reNode = outgroups.get(0);
		} else {
			reNode = EvolTreeOperator.getMostRecentCommonAncestor(outgroups.get(0), outgroups.get(1));
		}

		if (reNode != null) {
			root = EvolTreeOperator.rerootByNode(root, reNode);

			FileUtils.writeStringToFile(new File(outputNwkFile.getAbsolutePath() + ".rerooted.nwk"),
					treeCoder.code(root));
		} else {
			FileUtils.writeStringToFile(outputNwkFile, unRerootedTreeString);
		}

		return root;
	}

	/**
	 * 用其它的距离法来建树
	 * 
	 * @title reconstructTree
	 * @createdDate 2020-10-20 19:55
	 * @lastModifiedDate 2020-10-20 19:55
	 * @author "yudalang"
	 * @since 1.7
	 * 
	 * @param calculateDistanceMatrix
	 * @param outNames
	 * @return
	 * @return EvolNode
	 */
	EvolNode reconstructTree(float[][] calculateDistanceMatrix, String[] outNames) {
//		NodeEGPSv1 tree = new Upgma().tree(NJ4SarsCov2MultiThread.floatArray2double(calculateDistanceMatrix), outNames);
        NodeEGPSv1 tree = new NJ().tree(NJ4SarsCov2MultiThread.floatArray2double(calculateDistanceMatrix), outNames);
        return colneNode(tree);
	}

	private EvolNode colneNode(NodeEGPSv1 node) {
		EvolNode node1 = new DefaultPhyNode();

		if (node.getLeafName() != null) {
			node1.setName(node.getLeafName());
		}
		node1.setLength(node.getBranch().getLength());

		int childCount = node.getChildCount();

		if (childCount > 0) {
			for (int i = 0; i < childCount; i++) {
				EvolNode colneNode = colneNode(node.getChildAt(i));
				node1.addChild(colneNode);
			}
		}
		return node1;
	}

}
