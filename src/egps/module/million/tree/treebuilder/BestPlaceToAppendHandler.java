package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;

import module.evoltre.mutation.IMutation4Rec;

import java.util.ArrayList;
import java.util.List;

public class BestPlaceToAppendHandler {

	private final List<RecordOfMinMismatchPair>[] topElement2record;
	private final int size;

	// 最多存储 10000个
	private final int limit = 10000;

	@SuppressWarnings("unchecked")
	public BestPlaceToAppendHandler(int size) {

		if (size < 1) {
			throw new IllegalArgumentException("The size should greater than 0");
		}
		this.size = size;
		topElement2record = new List[size];
		for (int i = 0; i < size; i++) {
			topElement2record[i] = new ArrayList<>(limit);
		}
	}

	private RecordOfMinMismatchPair getRecord(boolean isMinMismatchBySplit, Node4AppendTree currNode,
			int mismatchToAdjIfBySplit, int minMismatch, List<IMutation4Rec> candidate2queryMuts,
			List<IMutation4Rec> splitMutationPartA, List<IMutation4Rec> splitMutationPartB,
			RecordOfMinMismatchPair record) {
		if (record == null) {
			return RecordOfMinMismatchPair.buildRecord(isMinMismatchBySplit, currNode, mismatchToAdjIfBySplit,
					minMismatch, candidate2queryMuts, splitMutationPartA, splitMutationPartB);
		} else {
			return record;
		}
	}

	/**
	 * 
	 * yudalang: 这里两个主要参数，minMismatch 和record。
	 * 如果已经直接生成了RecordOfMinMismatchPair对象，那么就不用再重新生成了。
	 * 使用这个方法来添加的话，不用先产生这个RecordOfMinMismatchPair对象。 先通过minMismatch来判断是否会被记录然后再生成对象.
	 * 
	 * mismatchToAdjIfBySplit 和 splitMutationPartA/splitMutationPartB 三个参数在
	 * isMinMismatchBySplit = false 的时候不起作用
	 * 
	 * @param minMismatch
	 */
	public void addOneEntery(int minMismatch, RecordOfMinMismatchPair record, 
			boolean isMinMismatchBySplit,
			Node4AppendTree currNode, int mismatchToAdjIfBySplit, List<IMutation4Rec> candidate2queryMuts,
			List<IMutation4Rec> splitMutationPartA, List<IMutation4Rec> splitMutationPartB) {

		for (int i = 0; i < size; i++) {
			List<RecordOfMinMismatchPair> currentList = topElement2record[i];

			if (currentList.isEmpty()) {
				RecordOfMinMismatchPair recordOfMinMismatchPair = getRecord(isMinMismatchBySplit, currNode,
						mismatchToAdjIfBySplit, minMismatch, candidate2queryMuts, splitMutationPartA,
						splitMutationPartB, record);
				currentList.add(recordOfMinMismatchPair);
				break;
			} else {
				RecordOfMinMismatchPair element = currentList.get(0);
				if (minMismatch == element.getMinMismatch()) {
					if (currentList.size() < limit) {
						RecordOfMinMismatchPair recordOfMinMismatchPair = getRecord(isMinMismatchBySplit, currNode,
								mismatchToAdjIfBySplit, minMismatch, candidate2queryMuts, splitMutationPartA,
								splitMutationPartB, record);
						currentList.add(recordOfMinMismatchPair);
					}
					break;
				} else if (minMismatch < element.getMinMismatch()) {

					List<RecordOfMinMismatchPair> last = topElement2record[size - 1];
					last.clear();
					RecordOfMinMismatchPair recordOfMinMismatchPair = getRecord(isMinMismatchBySplit, currNode,
							mismatchToAdjIfBySplit, minMismatch, candidate2queryMuts, splitMutationPartA,
							splitMutationPartB, record);
					last.add(recordOfMinMismatchPair);

					
					int next = size - 1;
					while (next > i) {
						topElement2record[next] = topElement2record[next - 1];
						next--;
					}

					topElement2record[i] = last;

					break;
				}
			}

		}
	}

	public List<RecordOfMinMismatchPair>[] getTopElement2record() {
		return topElement2record;
	}

	public void clear() {
		for (int i = 0; i < size; i++) {
			List<RecordOfMinMismatchPair> currentList = topElement2record[i];
			currentList.clear();
		}
	}

	public void addOneEntery(RecordOfMinMismatchPair recordOfMinMismatchPair) {
		addOneEntery(recordOfMinMismatchPair.getMinMismatch(), recordOfMinMismatchPair, false, null, limit, null, null, null);
		
	}
}
