package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.treeunit;

public class DefaultTreeNodeUnit implements TreeNodeUnit {

	@Override
	public TreeNodeUnit produceCluster(TreeNodeUnit b) {
		return new DefaultTreeNodeUnit();
	}

}
