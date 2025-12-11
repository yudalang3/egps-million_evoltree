package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit;

import java.util.Comparator;
import java.util.Objects;
import java.util.Optional;

public class SarsCov2TreeUnitPairComparatorSimple implements Comparator<TreeNodeUnitPair> {
	
	private static final int TIME_DIFF_TOLERATE = 30;

	public SarsCov2TreeUnitPairComparatorSimple() {
		
	}
	/**
	 * o1 better return 1; o2 better return -1; equal return 0
	 */
	@Override
	public int compare(TreeNodeUnitPair o1, TreeNodeUnitPair o2) {

		int compare;
		SarsCov2TreeUnitSimple firPairA = (SarsCov2TreeUnitSimple) o1.pairNodeUnitA;
		SarsCov2TreeUnitSimple firPairB = (SarsCov2TreeUnitSimple) o1.pairNodeUnitB;
		SarsCov2TreeUnitSimple secPairA = (SarsCov2TreeUnitSimple) o2.pairNodeUnitA;
		SarsCov2TreeUnitSimple secPairB = (SarsCov2TreeUnitSimple) o2.pairNodeUnitB;

		// 分别获取两对的时间差
		double firPairTime = Math.abs(firPairA.dayAfterEarliestDate - firPairB.dayAfterEarliestDate);
		double secPairTime = Math.abs(secPairA.dayAfterEarliestDate - secPairB.dayAfterEarliestDate);

		boolean firstMeetTimeTolerate = firPairTime < TIME_DIFF_TOLERATE;
		boolean secondMeetTimeTolerate = secPairTime < TIME_DIFF_TOLERATE;
		
		
		if (firstMeetTimeTolerate) {
			if (secondMeetTimeTolerate) {

			} else {
				//第一对的时间满足，第二对不满足时，第一对更好
				return 1;
			}
		} else {
			if (secondMeetTimeTolerate) {
				//第一对时间不满足，第二对时间满足时，第二对更好
				return -1;
			} else {

			}
		}
		
		// 当都满足时间的条件后。比较国家
		
		Optional<String> firPairString = getPairString(firPairA.country, firPairB.country);
		Optional<String> secPairString = getPairString(secPairA.country, secPairB.country);
		
		if (firPairString.isPresent()) {
			if (secPairString.isPresent()) {
				
			}else {
				//第一对的国家相同满足，第二对不相同时，第一对更好
				return 1;
			}
		}else {
			if (secPairString.isPresent()) {
				//第一对国家不相同，第二对国家相同时，第二对更好
				return -1;
			}else {
				
			}
		}
		
		//第一第二对时间都在容忍范围内，且国家都相同，再回看时间更近的
		if(firPairTime<secPairTime) {
			compare = 1;
		}else if(firPairTime>secPairTime){
			compare = -1;
		}
		
		compare=0;
		
//		// 先比较时间
//		// 值越小越好
//		if (secPairTime == firPairTime) {// 时间相等时比较国家
//			Optional<String> firPairString = getPairString(firPairA.country, firPairB.country);
//			Optional<String> secPairString = getPairString(secPairA.country, secPairB.country);
//			
//			if (firPairString.isPresent()) {
//				compare = 1;
//			}else {
//				if (secPairString.isPresent()) {
//					compare = -1;
//				}
//			}
//			compare = 0;
//		} else if (secPairTime < firPairTime) {
//			compare = -1;
//		} else {
//			compare = 1;
//		}
		return compare;
	}
	
	private Optional<String> getPairString(String a,String b) {
		String ret = null;
		if (Objects.nonNull(a) && Objects.nonNull(b)) {
			if (a.equals(b)) {
				ret = a;
			}
		}
		
		return Optional.ofNullable(ret);
	}

}
