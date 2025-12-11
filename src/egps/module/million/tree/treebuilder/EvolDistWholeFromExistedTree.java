package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;

import org.apache.commons.lang3.tuple.Pair;

/**
 * 
 *
 * @Description: Nothing change here, because default implementation is for
 *               this.
 *
 *
 * @createdDate Sep 27, 2020
 * @lastModifiedDate Sep 27, 2020
 * @author "yudalang"
 */
public class EvolDistWholeFromExistedTree extends AbstractEvolDistWholeRunnerSingleThread {

	@Override
	protected void calculateDistanceOfOutgroup2innerSequences() {
		final int refGenomeIndex = 0;
		/**
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
		for (Pair<String, Float> twoTuple : calDistanceFromOutgroups2ancestor) {
			float distOfOutgroup2ancestor = twoTuple.getRight();
			int outGroupIndex = seqName2seqNameIndexMap.get(twoTuple.getLeft());
			
			for (float[] fs : ret) {
				fs[outGroupIndex] = distOfOutgroup2ancestor + fs[refGenomeIndex];
			}
		}
		// Don't need to calculate the distance of outgroup1 to outgroup2.
	}

	@Override
	protected void initializeFinalDistanceMatrixVariable() {
		// new added sequences in the last!
		int[] js = sequenceNameIndexes[sequenceNameIndexes.length - 1];
		// This array not include reference genome!
		ret = new float[js.length - 1][seqName2seqNameIndexMap.size()];
	}

	/**
	 * 
	 * <p>
	 * Title: addDifference2matrix
	 * </p>
	 * <p>
	 * Description:
	 * 
	 * 
	 * 
	 * </p>
	 * 
	 * @param int1
	 * @param int2
	 * @param numOfDiff
	 * @see AbstractEvolDistWholeRunnerSingleThread#addDifference2matrix(int,
	 *      int, double)
	 *
	 */
	@Override
	protected void addDifference2matrix(int int1, int int2, double numOfDiff) {
		int numberOfAlreadyExistedOTUs = getNumberOfAlreadyExistedOTUs();
		
		int bigger, smaller;
		if (int1 > int2) {
			bigger = int1;smaller = int2;
		} else {
			bigger = int2;smaller = int1;
		}

		int indexOfRow = bigger - numberOfAlreadyExistedOTUs; 
		int indexOfCol = smaller;

		if (indexOfRow > -1) {
			ret[indexOfRow][indexOfCol] += numOfDiff;
		}

	}
}
