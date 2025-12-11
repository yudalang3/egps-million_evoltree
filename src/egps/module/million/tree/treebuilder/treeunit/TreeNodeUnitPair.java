package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit;

public final class TreeNodeUnitPair {

	int minI;
	int minJ;
	
	TreeNodeUnit pairNodeUnitA;
	TreeNodeUnit pairNodeUnitB;
	
	@Override
	public String toString() {
		return minI+"\t"+minJ + "\t"+pairNodeUnitA +"\t"+pairNodeUnitB+"\t";
	}


}
