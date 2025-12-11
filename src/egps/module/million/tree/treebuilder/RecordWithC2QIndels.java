package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;

import module.evoltre.mutation.IMutation4Rec;

import java.util.List;


public class RecordWithC2QIndels {

	RecordOfMinMismatchPair record;
	
	List<IMutation4Rec> c2qIndels;
	
	/**
	 * 一般情况下只有一个突变会是多重击中，
	 * 所以这个值一般是0或者1
	 */
	int indelMultipleHitCount;
	
	int getNumOfC2qIndels() {
		return c2qIndels.size();
	}
}
