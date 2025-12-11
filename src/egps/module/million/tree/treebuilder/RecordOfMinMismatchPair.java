package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;

import module.evoltre.mutation.IMutation4Rec;

import java.util.List;
import java.util.StringJoiner;


public class RecordOfMinMismatchPair {

	private Node4AppendTree targetNode;

	private boolean getMinMismatchBySplit;

	/**
	 * 不能简单认为 minMismatch = mutations.size() 因为有裂开的情况
	 */
	private int minMismatch;

	/**
	 * whatever getMinMismatchBySplit is,this is candidate to query mutations
	 */
	private List<IMutation4Rec> candidate2queryMutations;

	/**
	 * 如果Q早于Candidate 这个值为正值；否则是负值。
	 */
	long dateDiffCandidateMinusAdjustQuery;

	/**
	 * 如果是通过裂开多突变的分枝得到的最佳的目标节点，当遇到mismatch一样的情况时。选择最佳位点时需要校正时间
	 * 
	 * <pre>
	 *  
	 *                    |-------Candidate  mut1
	 *     |-------aa-----| 
	 *     |              |-----bb-----Q mut2
	 * ----|
	 *     |
	 *     |-----
	 *     
	 *     校正的Q相对于T的 mismatch应该是  #mut2 - #mut1
	 *     这个值有可能是负值，这表示应该加上更多的时间
	 * 
	 * </pre>
	 * 
	 * 如果不是通过通过裂开的话， splitMutationPartB 应该存储 P->Q的突变
	 */
	int mismatchToAdjIfBySplit;

	/**
	 * 裂开C的突变的时候，对应上面那幅图中的aa
	 */
	private List<IMutation4Rec> mutsAAorNone;
	/**
	 * 裂开C的突变的时候，对应上面那幅图中的bb。
	 * 否则就用来存储 P->Q的突变
	 */
	private List<IMutation4Rec> mutsBBOrP2Q;
	
	/**
	 * 在追加的过程中，我们有可能需要重新计算突变。
	 * 有没有模糊碱基的时候，我们可以复用 mutsAAorNone和mutsBBOrP2Q
	 */
	private boolean shouldReInferMutationsInAppending = false;
	
	RecordOfMinMismatchPair() {
		
	}

	public static RecordOfMinMismatchPair buildRecord(boolean isMinMismatchBySplit, Node4AppendTree currNode,
			int mismatchToAdjIfBySplit, int currMismatch, List<IMutation4Rec> candidate2queryMuts,
			List<IMutation4Rec> splitMutationPartA, List<IMutation4Rec> splitMutationPartB) {
		
		RecordOfMinMismatchPair recordOfMinMismatchPair = new RecordOfMinMismatchPair();
		recordOfMinMismatchPair.setGetMismatchToAdjIfBySplit(mismatchToAdjIfBySplit);
		recordOfMinMismatchPair.setMinMismatch(currMismatch);
		recordOfMinMismatchPair.setMutations(candidate2queryMuts);
		recordOfMinMismatchPair.setTargetNode(currNode);
		
		if (isMinMismatchBySplit) {
			if (splitMutationPartA == null || splitMutationPartB == null) {
				throw new IllegalArgumentException("When the best target is split, the splitMutationPartA and PartB can not be null.");
			}
		}
		recordOfMinMismatchPair.setGetMinMismatchBySplit(isMinMismatchBySplit);
		recordOfMinMismatchPair.setSplitMutationPartA(splitMutationPartA);
		recordOfMinMismatchPair.setSplitMutationPartB(splitMutationPartB);
		return recordOfMinMismatchPair;
	}

	@Override
	public String toString() {
		StringJoiner sbJoiner = new StringJoiner(" ");

		sbJoiner.add("ID:").add(targetNode.getCGBID());
		sbJoiner.add("Split:").add(String.valueOf(getMinMismatchBySplit));
		sbJoiner.add("MinMis:").add(String.valueOf(minMismatch));
		// 直接说明突变
		// sbJoiner.add("Muts:").add(candidate2queryMutations.toString());

		return sbJoiner.toString();
	}

	public Node4AppendTree getTargetNode() {
		return targetNode;
	}

	public void setTargetNode(Node4AppendTree targetNode) {
		this.targetNode = targetNode;
	}

	public boolean isGetMinMismatchBySplit() {
		return getMinMismatchBySplit;
	}

	public void setGetMinMismatchBySplit(boolean getMinMismatchBySplit) {
		this.getMinMismatchBySplit = getMinMismatchBySplit;
	}

	public int getMinMismatch() {
		return minMismatch;
	}

	public void setMinMismatch(int minMismatch) {
		this.minMismatch = minMismatch;
	}

	public List<IMutation4Rec> getMutations() {
		return candidate2queryMutations;
	}

	public void setMutations(List<IMutation4Rec> mutations) {
		this.candidate2queryMutations = mutations;
	}

	public long getDateDiffCandidateMinusAdjustQuery() {
		return dateDiffCandidateMinusAdjustQuery;
	}

	public void setDateDiffCandidateMinusAdjustQuery(long dateDiffCandidateMinusAdjustQuery) {
		this.dateDiffCandidateMinusAdjustQuery = dateDiffCandidateMinusAdjustQuery;
	}

	public void setGetMismatchToAdjIfBySplit(int mismatchToAdjIfBySplit) {
		this.mismatchToAdjIfBySplit = mismatchToAdjIfBySplit;

	}

	public void setSplitMutationPartA(List<IMutation4Rec> splitMutationPartA) {
		this.mutsAAorNone = splitMutationPartA;
	}

	public void setSplitMutationPartB(List<IMutation4Rec> splitMutationPartB) {
		this.mutsBBOrP2Q = splitMutationPartB;
	}

	public List<IMutation4Rec> getSplitMutationPartA() {
		return mutsAAorNone;
	}

	public List<IMutation4Rec> getSplitMutationPartB() {
		return mutsBBOrP2Q;
	}

	/**
	 * @return the shouldReInferMutationsInAppending
	 */
	public boolean isShouldReInferMutationsInAppending() {
		return shouldReInferMutationsInAppending;
	}

	/**
	 * @param shouldReInferMutationsInAppending the shouldReInferMutationsInAppending to set
	 */
	public void setShouldReInferMutationsInAppending(boolean shouldReInferMutationsInAppending) {
		this.shouldReInferMutationsInAppending = shouldReInferMutationsInAppending;
	}
	
}
