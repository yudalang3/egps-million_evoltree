package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import fasta.io.FastaReader;
import module.evoldist.operator.util.QuickDistUtil;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.tuple.Pair;
import phylo.msa.util.MsaCommonUtil;
import utils.EGPSUtil;


public class AbstractEvolDistWholeRunnerSingleThread {

	protected Map<String, Integer> seqName2seqNameIndexMap;

	String refGenomeID = "";
	String refGenome;

	protected float[][] ret;
	private List<String[]> inputSequences;
	protected int[][] sequenceNameIndexes;

	private List<File> inputFastaFiles;
	private File ancesterStateFile;
	private File fastaFileContainOutgroupFile;

	private int realStartPosition = 0;
	private int realEndPostion = 0;

	private int numOfAlreadyExistedOTUs;

	private Map<Float, String> ancestorStates;
	protected List<Pair<String, Float>> calDistanceFromOutgroups2ancestor;
	protected LinkedHashMap<String, String> outGroupSequenceMap;

//	DiversityCalculator diversityCalculator = new DiversityCalculator();

	protected AbstractEvolDistWholeRunnerSingleThread() {

	}

	public void setInputFastaFiles(List<File> files) {
		this.inputFastaFiles = new ArrayList<>();
		for (File file : files) {
			if ('.' == file.getName().charAt(0)) {
				System.err.println("Warning: ignore file: " + file.getAbsolutePath());

				continue;
			}
			this.inputFastaFiles.add(file);
		}

	}

	/**
	 * @param fastaFileContainOutgroupFile the {@link #fastaFileContainOutgroupFile}
	 *                                     to set
	 */
	public void setFastaFileContainOutgroupFile(File fastaFileContainOutgroupFile) {
		this.fastaFileContainOutgroupFile = fastaFileContainOutgroupFile;
	}

	public void setRefGenomeFile(File refGenomeFile) throws IOException {
		LinkedHashMap<String, String> readFastaDNASequence = FastaReader.readFastaDNASequence(refGenomeFile);
		Entry<String, String> next = readFastaDNASequence.entrySet().iterator().next();
		refGenome = next.getValue();
	}

	public void setRealCalculatePostion(int start, int end) {
		realStartPosition = start;
		realEndPostion = end;
	}

	public void setAncesterStateFile(File ancesterStateFile) {
		this.ancesterStateFile = ancesterStateFile;
	}

	public void setRefGenomeID(String refGenomeID) {
		this.refGenomeID = refGenomeID;
	}

	public float[][] calculateDistanceMatrix() throws IOException {

		// initialize and assign values
		int numOfFastaFiles = inputFastaFiles.size();
		sequenceNameIndexes = new int[numOfFastaFiles][];
		inputSequences = new ArrayList<String[]>(numOfFastaFiles);

		OneFileDistanceCalculator[] oneFileDistanceCalculators = new OneFileDistanceCalculator[numOfFastaFiles];

		// don't forget ancestor states

		// 赋值祖先状态
		ancestorStates = new HashMap<>(30000);
		int length = refGenome.length();
		for (int i = 0; i < length; i++) {
			ancestorStates.put(new Float(i + 1), String.valueOf(refGenome.charAt(i)));
		}
		if (ancesterStateFile != null) {
			List<String> lines = FileUtils.readLines(ancesterStateFile);
			for (String string : lines) {
				if (string.startsWith("#")) {
					continue;
				}
				String trimmedStr = string.trim();
				String[] split = trimmedStr.split("\t");
				Float position = new Float(split[0]);
				String stateString = split[1];
				ancestorStates.put(position, stateString);
			}
		}

		seqName2seqNameIndexMap = new HashMap<>(16384);
		seqName2seqNameIndexMap.put(refGenomeID, 0);
		int totalNamesIndex = 1;

		/**
		 * 如果这个距离被利用到从已有的树添加序列。那么要注意不能将外群放到最后，因为新加的序列与原有的序列要满足这样的排列关系：
		 * 
		 * <pre>
		 * The matrix ret is :
		 *                  oldSeq1     oldSeq2    ...   outGroup1 outGroup2 (If out-group existed.)  newAddedSeq1  newAddedSeq2  newAddedSeq3  ... 
		 * newAddedSeq1
		 * newAddedSeq2
		 * newAddedSeq3
		 *     .
		 *     .
		 *     .
		 * </pre>
		 */
		for (int k = 0; k < numOfFastaFiles; k++) {
			File inputFastaFile = inputFastaFiles.get(k);

			LinkedHashMap<String, String> readFastaDNASequence = FastaReader.readFastaDNASequence(inputFastaFile,
					true, true);

			if (TreeBuilder4SarsCov2.debug) {
				System.out.println(
						"Num of sequences in file " + inputFastaFile.getName() + " is " + readFastaDNASequence.size());
			}

			List<Integer> nameIndexesOfOneFile = new LinkedList<>();
			List<String> listOfSequencesOfOneFile = new LinkedList<>();
			int index = 0;
			for (Entry<String, String> entry : readFastaDNASequence.entrySet()) {
				String key = entry.getKey();
				String value = entry.getValue();

				if (refGenomeID.equals(key)) {
					listOfSequencesOfOneFile.add(0, value);
					nameIndexesOfOneFile.add(0, 0);
					index++;
				} else {
					if (seqName2seqNameIndexMap.get(key) == null) {

						seqName2seqNameIndexMap.put(key, totalNamesIndex);

						nameIndexesOfOneFile.add(totalNamesIndex);

						listOfSequencesOfOneFile.add(value);
						index++;
						totalNamesIndex++;
					} else {
						System.err.println("Warning: " + key + " in file: " + inputFastaFile.getName() + " repeated!");
					}

				}
			}

			if (TreeBuilder4SarsCov2.debug && (k == numOfFastaFiles - 1)) {
				System.out.println("Num of non-repeated sequences in file (include refGenome) is: " + index);
			}
			readFastaDNASequence = null;

			int[] ints = nameIndexesOfOneFile.stream().mapToInt(Integer::valueOf).toArray();
			sequenceNameIndexes[k] = ints;
			String[] array = listOfSequencesOfOneFile.toArray(new String[2]);
			inputSequences.add(array);

			OneFileDistanceCalculator oneFileDistanceCalculator = new OneFileDistanceCalculator(array, ints);
			oneFileDistanceCalculator.setRunner(this);
			oneFileDistanceCalculators[k] = oneFileDistanceCalculator;

			// numOfFastaFiles - 1 is the last fasta file index.
			if (numOfFastaFiles == 1) {
				// only one file
				if (fastaFileContainOutgroupFile != null) {
					calDistanceFromOutgroups2ancestor(fastaFileContainOutgroupFile);

					for (Pair<String, Float> twoTuple : calDistanceFromOutgroups2ancestor) {
						seqName2seqNameIndexMap.put(twoTuple.getLeft(), totalNamesIndex++);
					}

				}
			} else {
				if (k == (numOfFastaFiles - 2)) {
					// 倒数第二个文件后面初始化这个
					if (fastaFileContainOutgroupFile != null) {
						calDistanceFromOutgroups2ancestor(fastaFileContainOutgroupFile);

						for (Pair<String, Float> twoTuple : calDistanceFromOutgroups2ancestor) {
							seqName2seqNameIndexMap.put(twoTuple.getLeft(), totalNamesIndex++);
						}
					}
				}
			}

		}

		System.out.println("Finishing loading fasta files!");

		// Calculate this value for the calculate distance of existed tree.
		// remove ref.
		int numOfNewAddedOTUs = sequenceNameIndexes[sequenceNameIndexes.length - 1].length - 1;
		numOfAlreadyExistedOTUs = seqName2seqNameIndexMap.size() - numOfNewAddedOTUs;

		initializeFinalDistanceMatrixVariable();

		if (TreeBuilder4SarsCov2.debug) {
			System.out.println(EGPSUtil.getAlreadyUsedJVMMemory() + " MB memory space has already use!");
		}

		// real running the process

		// 在计算每个位点时临时存储不同的状态对应的序列名称索引
		Map<SNPAndInsertionDiffTypeAspect, Collection<Integer>> onePositionDIffType2NameIndexMap = new HashMap<>();

		// 遍历处理每个位点
		for (int i = realStartPosition; i <= realEndPostion; i++) {

			// 遍历处理每个fasta文件的当前位点
			for (OneFileDistanceCalculator oneFileDistanceCalculator : oneFileDistanceCalculators) {
				oneFileDistanceCalculator.setRealCalculatedPosition(i);
				oneFileDistanceCalculator.processOneFile(onePositionDIffType2NameIndexMap);
			}

			// 将位点状态对应序列名的列表转换成距离矩阵
			listOfTotalMap2DistMatrix(onePositionDIffType2NameIndexMap, i);
			onePositionDIffType2NameIndexMap.clear();

		}

		/**
		 * 将病毒序列到外群的距离赋值上去：
		 */
		if (fastaFileContainOutgroupFile != null) {
			calculateDistanceOfOutgroup2innerSequences();
		}

		return ret;
	}

	/**
	 * 当外群文件存在时，计算外群当内群各个OUT之间的距离。
	 * 
	 * @title calculateDistanceOfOutgroup2innerSequences
	 * @createdDate 2020-11-05 09:20
	 * @lastModifiedDate 2020-11-05 09:20
	 * @author yudalang
	 * @since 1.7
	 * 
	 * @return void
	 */
	protected void calculateDistanceOfOutgroup2innerSequences() {
		final int refGenomeIndex = 0;
		final int totalNum = ret.length;
		for (Pair<String, Float> twoTuple : calDistanceFromOutgroups2ancestor) {
			float dist = twoTuple.getRight();
			int outGroupIndex = seqName2seqNameIndexMap.get(twoTuple.getLeft());

			for (int i = 0; i < outGroupIndex; i++) {
				if (outGroupIndex > i) {
					ret[outGroupIndex][i] = dist + ret[i][refGenomeIndex];
				} else {
					ret[i][outGroupIndex] = dist + ret[i][refGenomeIndex];
				}
			}

			for (int i = outGroupIndex + 1; i < totalNum; i++) {
				if (outGroupIndex > i) {
					ret[outGroupIndex][i] = dist + ret[i][refGenomeIndex];
				} else {
					ret[i][outGroupIndex] = dist + ret[i][refGenomeIndex];
				}
			}
		}

		int size = calDistanceFromOutgroups2ancestor.size();
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < i; j++) {
                Pair<String, Float> a = calDistanceFromOutgroups2ancestor.get(i);
                Pair<String, Float> b = calDistanceFromOutgroups2ancestor.get(j);
				int index1 = seqName2seqNameIndexMap.get(a.getLeft());
				String seq1 = outGroupSequenceMap.get(a.getLeft());
				int index2 = seqName2seqNameIndexMap.get(b.getLeft());
				String seq2 = outGroupSequenceMap.get(b.getLeft());

				float d = (float) MsaCommonUtil.MSAUtil.getPairwiseSequenceDistance(seq1, seq2);
				// System.err.println(index1 + "\t" + index2 + "\t" + d);

				if (index1 > index2) {
					ret[index1][index2] = d;
				} else {
					ret[index2][index1] = d;
				}
			}
		}
	}

	private void calDistanceFromOutgroups2ancestor(File file) throws IOException {
		outGroupSequenceMap = FastaReader.readFastaDNASequence(file);

		String refGenome = outGroupSequenceMap.get(refGenomeID);
		outGroupSequenceMap.remove(refGenomeID);

		List<String> list = new ArrayList<String>(outGroupSequenceMap.values());
		List<String> seqNames = new ArrayList<String>(outGroupSequenceMap.keySet());

		int size = list.size();
		List<Pair<String, Float>> ret = new ArrayList<>();
		for (int i = 0; i < size; i++) {
			double pairwiseSequenceDistance = MsaCommonUtil.MSAUtil.getPairwiseSequenceDistance(list.get(i), refGenome);
			String string = seqNames.get(i);
			ret.add(Pair.of(string, new Float(pairwiseSequenceDistance)));
		}

		calDistanceFromOutgroups2ancestor = ret;

	}

	/**
	 * 
	 * 初始化的过程会计算祖先序列与外群之间的距离，同时再把外群的数量统计出来。
	 * 
	 * @title initializeFinalDistanceMatrixVariable
	 * @createdDate 2020-10-20 09:37
	 * @lastModifiedDate 2020-10-20 09:37
	 * @author "yudalang"
	 * @since 1.7
	 * 
	 * @return void
	 * @throws IOException
	 */
	protected void initializeFinalDistanceMatrixVariable() throws IOException {

		int numOfTotalOTU = seqName2seqNameIndexMap.size();
		ret = new float[numOfTotalOTU][];
		for (int i = 0; i < numOfTotalOTU; i++) {
			ret[i] = new float[i + 1];
		}
	}

	Set<Character> atcgs = new HashSet<>(Arrays.asList('A', 'T', 'C', 'G', 'a', 't', 'c', 'g'));

	private void listOfTotalMap2DistMatrix(Map<SNPAndInsertionDiffTypeAspect, Collection<Integer>> mapOfPos, int key)
			throws IOException {

		// 没有多态性
		int size = mapOfPos.size();
		if (size < 2) {
			return;
		}

		// don't forget the ref genome
		// Since we use the HashSet, so we don't need to add

		List<SNPAndInsertionDiffTypeAspect> iterater = new ArrayList<>(mapOfPos.keySet());
		List<Collection<Integer>> valuesOfCollections = new ArrayList<>(mapOfPos.values());

		if (TreeBuilder4SarsCov2.debug) {
			System.out.println("Position: " + key + " num. of DiffTypeWithTwoAspects is " + iterater.size() + "  ");
//			for (Entry<SNPAndInsertionDiffTypeAspect, Collection<Integer>> entry : mapOfPos.entrySet()) {
//				System.out.println(entry);
//			}
		}

		coreProcess(iterater, key, valuesOfCollections);

		return;
	}

	private void coreProcess(List<SNPAndInsertionDiffTypeAspect> iterater, int key,
			List<Collection<Integer>> valuesOfCollections) {

		int sizeOfType = iterater.size();
		// 保险性判定
		if (sizeOfType == 1) {
			return;
		}

		for (int i = 0; i < sizeOfType; i++) {
			for (int j = i + 1; j < sizeOfType; j++) {
				SNPAndInsertionDiffTypeAspect diffType1 = iterater.get(i);
				SNPAndInsertionDiffTypeAspect diffType2 = iterater.get(j);

				// System.out.println("i= "+i+" j= "+j+" total= "+sizeOfType);

				Collection<Integer> set1 = valuesOfCollections.get(i);
				Collection<Integer> set2 = valuesOfCollections.get(j);
				double cc = getDistanceOfTwoDiffType(diffType1, diffType2);
				if (cc > 0) {
					for (Integer int1 : set1) {
						for (Integer int2 : set2) {

							if (int1 == int2) {
								System.err.println("Has the same index in position: " + key);
								// System.out.println(mapOfPos);
							} else {
								addDifference2matrix(int1, int2, cc);
							}
						}
					}
				}
			}
		}

	}

    public static double getDistanceOfTwoDiffType(SNPAndInsertionDiffTypeAspect diffType1,
                                                  SNPAndInsertionDiffTypeAspect diffType2) {

        double ret = 0;

        // insitu place
        if (diffType1.isDeletion()) {
            if (diffType2.isDeletion()) {
                if (diffType1.isIfFirstDeletion() && diffType2.isIfFirstDeletion()) {
                    if (diffType1.getLengthOfDeletionRefer2RefSeq() == diffType2.getLengthOfDeletionRefer2RefSeq()) {
                        ret++;
                    }
                }
            } else {
                if (diffType1.isIfFirstDeletion()) {
                    ret++;
                }
            }
        } else {
            if (diffType2.isDeletion()) {
                if (diffType2.isIfFirstDeletion()) {
                    ret++;
                }
            } else {
                ret += QuickDistUtil.getTwoSNPCharDifferenceWithAmbiguousBaseAccording2IntArray(
                        diffType1.getInsituSite(), diffType2.getInsituSite());
            }
        }

        // insetion String
        boolean type1IsNotRight = diffType1.getInsertionContent().isEmpty() && diffType1.isDeletion();
        boolean type2IsNotRight = diffType2.getInsertionContent().isEmpty() && diffType2.isDeletion();

        if (type1IsNotRight || type2IsNotRight) {
            // 其中之一有问题都不行
        } else {
            //			if (!diffType1.getInsertionContent().equals(diffType2.getInsertionContent())) {
            //				ret ++;
            //			}
            // 考虑模糊碱基
            if (!QuickDistUtil.judgeTwoAllelesIdentities(diffType1.getInsertionContent(),
                    diffType2.getInsertionContent())) {
                ret++;
            }
        }
        return ret;
    }

	protected void addDifference2matrix(int int1, int int2, double numOfDiff) {
		if (int1 > int2) {
			ret[int1][int2] += numOfDiff;
		} else {
			ret[int2][int1] += numOfDiff;
		}
	}

	public String[] getOUTNames() {
		List<Map.Entry<String, Integer>> list = new LinkedList<>(seqName2seqNameIndexMap.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {
			@Override
			public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
				return (o1.getValue()).compareTo(o2.getValue());
			}
		});

		int size = list.size();
		String[] ret = new String[size];

		for (int i = 0; i < size; i++) {
			ret[i] = list.get(i).getKey();
		}

		return ret;
	}

	Collection<Integer> newInstanceOfDiffType2CollectionMapVaules() {
		Collection<Integer> ret = new HashSet<Integer>(1024);
		return ret;
	}

	Collection<Integer> newInstanceOfDiffType2CollectionMapVaules(Collection<Integer> input) {
		Collection<Integer> ret = new HashSet<Integer>(input);
		return ret;
	}

	public int getNumberOfAlreadyExistedOTUs() {
		return numOfAlreadyExistedOTUs;
	}

	public Map<String, Integer> getSeqName2seqNameIndexMap() {
		return seqName2seqNameIndexMap;
	}

	public String[] getLastAlignment() {
		return inputSequences.get(inputSequences.size() - 1);
	}
}
