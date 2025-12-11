package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;


import java.util.Objects;

import static phylo.msa.util.EvolutionaryProperties.GAP_CHAR;

public class SNPAndInsertionDiffTypeAspect {

	private final char insituSite;
	private final String insertionContent;
	private final int lengthOfDeletionRefer2RefSeq;
	/**
	 * 当 insituSite 是 ‘-’ 的时候，这个参数就有用了。
	 */
	private final boolean ifFirstDeletion;

	/**
	 * 
	 * For SNP
	 *  
	 * @title SNPAndInsertionDiffTypeAspect
	 * @createdDate 2020-10-25 11:23
	 * @lastModifiedDate 2020-10-25 11:23
	 * @author yudalang
	 * @since 1.7
	 *   
	 * @param insituSite
	 */
	public SNPAndInsertionDiffTypeAspect(char insituSite){
		this(insituSite, "", 0, false);
	}
	
	/**
	 * 
	 * For not first deletion.
	 *  
	 * @title SNPAndInsertionDiffTypeAspect
	 * @createdDate 2020-10-25 11:23
	 * @lastModifiedDate 2020-10-25 11:23
	 * @author yudalang
	 * @since 1.7
	 *   
	 * @param insituSite
	 * @param insertionContent
	 */
	public SNPAndInsertionDiffTypeAspect(char insituSite,String insertionContent){
		this(insituSite, insertionContent, 0, false);
	}
	
	/**
	 * 
	 * for first deletion.
	 *  
	 * @title SNPAndInsertionDiffTypeAspect
	 * @createdDate 2020-10-25 11:23
	 * @lastModifiedDate 2020-10-25 11:23
	 * @author yudalang
	 * @since 1.7
	 *   
	 * @param insituSite
	 * @param insertionContent
	 * @param lengthOfDeletionRefer2RefSeq
	 * @param ifFirstDeletion
	 */
	public SNPAndInsertionDiffTypeAspect(char insituSite, 
			String insertionContent, 
			int lengthOfDeletionRefer2RefSeq,
			boolean ifFirstDeletion) {
		
		this.insituSite = insituSite;
		this.insertionContent = insertionContent;
		this.lengthOfDeletionRefer2RefSeq = lengthOfDeletionRefer2RefSeq;
		this.ifFirstDeletion = ifFirstDeletion;
	}

	public char getInsituSite() {
		return insituSite;
	}

	public String getInsertionContent() {
		return insertionContent;
	}

	public int getLengthOfDeletionRefer2RefSeq() {
		return lengthOfDeletionRefer2RefSeq;
	}

	public boolean isIfFirstDeletion() {
		return ifFirstDeletion;
	}

	public boolean isDeletion() {
		return this.insituSite == GAP_CHAR;
	}

	@Override
	public int hashCode() {
		return Objects.hash(ifFirstDeletion, insertionContent, insituSite);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (!(obj instanceof SNPAndInsertionDiffTypeAspect)) {
			return false;
		}
		SNPAndInsertionDiffTypeAspect other = (SNPAndInsertionDiffTypeAspect) obj;
		return ifFirstDeletion == other.ifFirstDeletion && Objects.equals(insertionContent, other.insertionContent)
				&& insituSite == other.insituSite;
	}

	@Override
	public String toString() {
		return "SNPAndInsertionDiffTypeAspect [insituSite=" + insituSite + ", insertionContent=" + insertionContent
				+ ", lengthOfDeletionRefer2RefSeq=" + lengthOfDeletionRefer2RefSeq + ", ifFirstDeletion="
				+ ifFirstDeletion + "]";
	}
	
	
	
}
