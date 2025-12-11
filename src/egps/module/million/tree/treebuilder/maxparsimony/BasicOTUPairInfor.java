package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.maxparsimony;

import module.evoltre.mutation.IMutation4Rec;

import java.util.List;


public class BasicOTUPairInfor {
	
	BasicOTUInfor seq1Infor;
	
	BasicOTUInfor seq2Infor;
	
	List<IMutation4Rec> commonMutations;
	int numOfCommonMutations;
	
	
	
	/**
	 * @return the {@link #seq1Infor}
	 */
	public BasicOTUInfor getSeq1Infor() {
		return seq1Infor;
	}
	/**
	 * @param seq1Infor the {@link #seq1Infor} to set
	 */
	public void setSeq1Infor(BasicOTUInfor seq1Infor) {
		this.seq1Infor = seq1Infor;
	}
	/**
	 * @return the {@link #seq2Infor}
	 */
	public BasicOTUInfor getSeq2Infor() {
		return seq2Infor;
	}
	/**
	 * @param seq2Infor the {@link #seq2Infor} to set
	 */
	public void setSeq2Infor(BasicOTUInfor seq2Infor) {
		this.seq2Infor = seq2Infor;
	}
	/**
	 * @return the {@link #commonMutations}
	 */
	public List<IMutation4Rec> getCommonMutations() {
		return commonMutations;
	}
	/**
	 * @param commonMutations the {@link #commonMutations} to set
	 */
	public void setCommonMutations(List<IMutation4Rec> commonMutations) {
		this.commonMutations = commonMutations;
	}
	/**
	 * @return the {@link #numOfCommonMutations}
	 */
	public int getNumOfCommonMutations() {
		return numOfCommonMutations;
	}
	/**
	 * @param numOfCommonMutations the {@link #numOfCommonMutations} to set
	 */
	public void setNumOfCommonMutations(int numOfCommonMutations) {
		this.numOfCommonMutations = numOfCommonMutations;
	}
	
	
	

}
