package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.multithread;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.TreeSet;

import egps.module.ncov.analysis.egps.module.million.tree.treebuilder.OTUState.OneLeafWithAllStates;
import module.parsimonytre.algo.StateAfterMutation;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.ml.distance.EuclideanDistance;

import utils.EGPSObjectCounter;

public class Function_ClusterSequences {

	

	public static Pair<List<String>, List<String>> getTwoGroupOfSequences4appending(List<String> inputSeqs)
			throws IOException {

		Comparator<Entry<Byte, MutableInt>> comparator = new Comparator<Entry<Byte, MutableInt>>() {
			@Override
			public int compare(Entry<Byte, MutableInt> o1, Entry<Byte, MutableInt> o2) {
				return o2.getValue().compareTo(o1.getValue());
			}
		};

		List<String> outputStrings = new LinkedList<>();

		File finalLeafStatesFile = new File(
				"/Users/yudalang/hugeDataRepo/sarsCov2fastaData/indexData/allLeafStates1102_1129.txt");
		List<String> asList = readOneFileContent();

		Pair<List<String>, List<String>> twoDistinctPairOfSeqs = Pair.of(new ArrayList<>(), new ArrayList<>());

        Map<String, OneLeafWithAllStates> finalLeafStatesMapOfInterest = new  HashMap<>();
        //TODO 需要填充

		double[][] data = getAllSites(finalLeafStatesMapOfInterest, asList);

		int exceptDividedSize = asList.size() / 2;
		int tolerate = (int) (0.1 * exceptDividedSize);
		if (tolerate < 1) {
			tolerate = 1;
		}

		int k4cluster = 2;
		KmeansPlusPlus kmeansPlusPlus = new KmeansPlusPlus(k4cluster, 200, new EuclideanDistance());
		EGPSObjectCounter<Byte> stringCounter = new EGPSObjectCounter<>();

		boolean suitableResult = false;
		while (!suitableResult) {
			kmeansPlusPlus.setK(k4cluster);
			kmeansPlusPlus.setMatrix(data);
			byte[] doClustering = kmeansPlusPlus.doClustering();

			stringCounter.clear();
			for (byte b : doClustering) {
				stringCounter.addOneEntry(b);
			}
			Map<Byte, MutableInt> counterMap = stringCounter.getCounterMap();
			List<Entry<Byte, MutableInt>> list = new ArrayList<>(counterMap.entrySet());

			Collections.sort(list, comparator);

			System.out.println("K is " + k4cluster);
			System.out.println("print cluster results:\n" + list.toString());

			k4cluster += 2;

			Entry<Byte, MutableInt> biggestEntry = list.get(0);
			if (k4cluster > 10 || Math.abs(exceptDividedSize - biggestEntry.getValue().intValue()) < tolerate) {
				byte chosenClass = biggestEntry.getKey();
				System.out.println("Chosen class is " + chosenClass + "\t" + biggestEntry.getValue());

				int index = 0;
				for (byte entry : doClustering) {
					String string = asList.get(index);
					if (entry == chosenClass) {
						twoDistinctPairOfSeqs.getLeft().add(string);
					} else {
						twoDistinctPairOfSeqs.getRight().add(string);
					}
					index++;
				}
				break;
			}
		}

		System.out.println(twoDistinctPairOfSeqs.getLeft().size() + "\t" + twoDistinctPairOfSeqs.getRight().size());
		System.out.println(twoDistinctPairOfSeqs.getLeft());
		System.out.println(twoDistinctPairOfSeqs.getRight());

		FileUtils.writeLines(new File("/Users/yudalang/Desktop/11.csv"), outputStrings);

		return twoDistinctPairOfSeqs;
	}

	private static List<String> readOneFileContent() {
		List<String> readLines = null;
		try {
			readLines = FileUtils.readLines(new File("/Users/yudalang/Desktop/11.txt"));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return readLines;
	}

	private static double[][] getAllSites(Map<String, OneLeafWithAllStates> finalLeafStatesMapOfInterest,
			List<String> sequences) {
		// 首先得到所有位置, TreeSet自动会对位置进行排序
		TreeSet<Integer> allSites = new TreeSet<>();
		for (Entry<String, OneLeafWithAllStates> entry : finalLeafStatesMapOfInterest.entrySet()) {
			List<StateAfterMutation> allStates = entry.getValue().getAllStates();
			for (StateAfterMutation state : allStates) {
				allSites.add(state.getPosition());
			}
		}

//		System.out.println("All sites are:");
//		System.out.println(allSites);
		int totalSite = allSites.size();
		int numOfSeqs = sequences.size();
		double[][] ret = new double[numOfSeqs][];

		// 产生一个基因组位置到当前位置的map
		int index = 0;
		Map<Integer, Integer> currentPos2indexMap = new HashMap<>(numOfSeqs);
		for (Integer ds : allSites) {
			currentPos2indexMap.put(ds, index);
			index++;
		}
		// 产生一个序列名字到以 0,1为代表的 haplotype 的map
		Map<String, double[]> seq2haptypeMap = new HashMap<>(numOfSeqs);
		for (Entry<String, OneLeafWithAllStates> entry : finalLeafStatesMapOfInterest.entrySet()) {

			double[] ds = new double[totalSite];
			List<StateAfterMutation> allStates = entry.getValue().getAllStates();
			for (StateAfterMutation state : allStates) {
				Integer integer = currentPos2indexMap.get(state.getPosition());
				ds[integer] = state.getInsutePlaceState();
			}

			seq2haptypeMap.put(entry.getKey(), ds);
		}

		// 最后根据序列名字将这个double数组生成
		index = 0;
		for (String name : sequences) {
			double[] ds = seq2haptypeMap.get(name);

			Objects.requireNonNull(ds, name + " can not get haplotype.");
			ret[index] = ds;
			index++;
		}

		return ret;
	}

	public static void main(String[] args) throws IOException {
		getTwoGroupOfSequences4appending(null);
	}

}
