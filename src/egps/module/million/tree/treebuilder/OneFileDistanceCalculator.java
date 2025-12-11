package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;


import phylo.msa.util.MsaCommonUtil;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import static phylo.msa.util.EvolutionaryProperties.GAP_CHAR;


public class OneFileDistanceCalculator {

	AbstractEvolDistWholeRunnerSingleThread runner;

	int realCalculatedPosition = 0;

	final String[] inputFastaCharArr;
	final int[] sequenceNameIndexes;

	int zeroBasedStartPos;
	int currentPositionWithoutGap = 1;

	public OneFileDistanceCalculator(String[] inputFastaCharArr, int[] sequenceNameIndexes) {
		super();
		this.inputFastaCharArr = inputFastaCharArr;
		this.sequenceNameIndexes = sequenceNameIndexes;

		initializeProperties();
	}

	private void initializeProperties() {
		String referenceGenomeSeq = inputFastaCharArr[0];
		int currPos = 0;
		while (referenceGenomeSeq.charAt(currPos) == GAP_CHAR) {
			currPos++;
		}
		zeroBasedStartPos = currPos;
	}

	public void setRunner(AbstractEvolDistWholeRunnerSingleThread runner) {
		this.runner = runner;
	}

	public void setRealCalculatedPosition(int realCalculatedPosition) {
		this.realCalculatedPosition = realCalculatedPosition;
	}

	public void processOneFile(Map<SNPAndInsertionDiffTypeAspect, Collection<Integer>> onePositionDIffType2NameIndexMap)
			throws IOException {

		String referenceGenomeSeq = inputFastaCharArr[0];
		int numOfTotalSites = referenceGenomeSeq.length();

		// get one position location
		int nextPos = 0;

		while (true) {
			
			nextPos = MsaCommonUtil.MSAUtil.getNextRefGenomePositionInAlignment(numOfTotalSites, referenceGenomeSeq, zeroBasedStartPos);
			
			if (nextPos == -1) {
				break;
			}
			
			if (currentPositionWithoutGap == realCalculatedPosition) {
				
//				if (TreeBuilder4SarsCov2.debug) {
//					System.out.println("Now, we start to deal with postion: " + currentPositionWithoutGap + "; currPos: "
//							+ zeroBasedStartPos + "\tnextPos: " + nextPos);
//				}
				
				processOnePositionAccording2RefGenome(referenceGenomeSeq, inputFastaCharArr, zeroBasedStartPos, nextPos,
						onePositionDIffType2NameIndexMap, sequenceNameIndexes);

				zeroBasedStartPos = nextPos;
				currentPositionWithoutGap++;
				break;
			}
			
			zeroBasedStartPos = nextPos;
			currentPositionWithoutGap++;
		}

	}

	private void processOnePositionAccording2RefGenome(String alignedRefGenome, String[] inputFastaCharArr,
			int inclusiveStartPosOnThisSite, int exclusiveEndPosOnThisSite,
			Map<SNPAndInsertionDiffTypeAspect, Collection<Integer>> onePositionDIffType2NameIndexMap,
			int[] sequenceNameIndexes) {

		int numOfSequence = inputFastaCharArr.length;

		for (int j = 0; j < numOfSequence; j++) {
			String currentGenomeSeq = inputFastaCharArr[j];
			char charAtHere = currentGenomeSeq.charAt(inclusiveStartPosOnThisSite);

			String insertionString = currentGenomeSeq.substring(inclusiveStartPosOnThisSite + 1,
					exclusiveEndPosOnThisSite);
			
			if (insertionString.indexOf( (int) GAP_CHAR) != -1) {
				insertionString = insertionString.replaceAll(String.valueOf(GAP_CHAR), "");
			}

			SNPAndInsertionDiffTypeAspect mutation = null;
			if (charAtHere == GAP_CHAR) {
				// Deletion
				if (ifFirstOccurredDeletion(alignedRefGenome, inclusiveStartPosOnThisSite, currentGenomeSeq)) {
					// first deletion
					// need to get the length of deletion ,only with refgenome

					int deletionLength = getDeletionLength(alignedRefGenome, currentGenomeSeq, inclusiveStartPosOnThisSite);
					if (deletionLength == 0) {
						// not first deletion
//						System.err.println("Has seq tail missing:");
//						System.err.println(currentGenomeSeq.substring(0, 1000));
						mutation = new SNPAndInsertionDiffTypeAspect(charAtHere,
								insertionString);
					}else {
						mutation = new SNPAndInsertionDiffTypeAspect(charAtHere, insertionString,
								deletionLength, true);
					}
				} else {
					// not first deletion
					mutation = new SNPAndInsertionDiffTypeAspect(charAtHere,
							insertionString);
				}
			} else {
				mutation = new SNPAndInsertionDiffTypeAspect(charAtHere, insertionString);

			}

			addEntery2HashMap(onePositionDIffType2NameIndexMap, mutation, sequenceNameIndexes[j]);

		}

	}

	private int getDeletionLength(String alignedRefGenome, String currentGenomeSeq, int inclusiveStartPosOnThisSite) {
		int ret = 0;
		int length = currentGenomeSeq.length();
		int endPos = MsaCommonUtil.MSAUtil.getNextRefGenomePositionInAlignment(length,
				currentGenomeSeq, inclusiveStartPosOnThisSite);

		
		if (endPos== -1) {
			endPos = length;
		}
		String substring = alignedRefGenome.substring(inclusiveStartPosOnThisSite, endPos);
		
		
		int ll = substring.length();

		for (int i = 0; i < ll; i++) {
			char charAt = substring.charAt(i);
			if (charAt != GAP_CHAR) {
				ret++;
			}
		}
		return ret;
	}

	/**
	 * 
	 * If all char is gap, this is not the first occurred mutation.
	 * <pre>
	 * So in case:
	 * Ref    ATCGACAGCTAGCTAGCATGCATCGATCGATCGAT...
	 * Seq1   -----------GCTAGCATGCATCGATCGATCGAT...
	 *                        | This is real process postion.
	 * </pre>
	 *  So it can handle this situation.
	 *  
	 * @title ifFirstOccurredDeletion
	 * @createdDate 2020-11-06 11:23
	 * @lastModifiedDate 2020-11-06 11:23
	 * @author yudalang
	 * @since 1.7
	 *   
	 * @param alignedRefGenome
	 * @param inclusiveStartPosOnThisSite
	 * @param currentGenomeSeq
	 * @return boolean
	 */
	private boolean ifFirstOccurredDeletion(String alignedRefGenome, int inclusiveStartPosOnThisSite,
			String currentGenomeSeq) {
		int previousRefGenomePositionInAlignment = MsaCommonUtil.MSAUtil.getPreviousRefGenomePositionInAlignment(alignedRefGenome,
				inclusiveStartPosOnThisSite);
		String substring = currentGenomeSeq.substring(previousRefGenomePositionInAlignment,
				inclusiveStartPosOnThisSite);

		int length = substring.length();
		for (int i = 0; i < length; i++) {
			char charAt = substring.charAt(i);
			if (charAt != GAP_CHAR) {
				return true;
			}
		}

		// if all char is gap, this is not the first occurred mutation.
		return false;
	}

	private void addEntery2HashMap(Map<SNPAndInsertionDiffTypeAspect, Collection<Integer>> snpHashMap,
			SNPAndInsertionDiffTypeAspect diffType, int seqNameIndex) {
		Collection<Integer> set = snpHashMap.get(diffType);
		if (set == null) {
			set = runner.newInstanceOfDiffType2CollectionMapVaules();
			set.add(seqNameIndex);
			snpHashMap.put(diffType, set);
		} else {
			set.add(seqNameIndex);
		}
	}

}
