package egps.module.ncov.analysis.egps.module.million.tree.treebuilder;

import module.evoltre.mutation.IMutation4Rec;
import module.evolview.model.tree.GraphicsNode;
import module.evolview.model.tree.NodeWithCGBID;
import module.parsimonytre.algo.Node4SankoffAlgo;
import module.parsimonytre.algo.StateAfterMutation;

import java.util.ArrayList;
import java.util.List;


public class Node4AppendTree extends NodeWithCGBID {

	private List<IMutation4Rec> cumulatedMutations;

	private List<StateAfterMutation> listOfStateAfterMutations;

	public Node4AppendTree() {
		this.cumulatedMutations = new ArrayList<>();
	}

	public void setCumulatedMutations(List<IMutation4Rec> cumulatedMutations) {
		this.cumulatedMutations = cumulatedMutations;
	}

	public List<IMutation4Rec> getCumulatedMutations() {
		return cumulatedMutations;
	}

	public void setListOfStateAfterMutations(List<StateAfterMutation> listOfStateAfterMutations) {
		this.listOfStateAfterMutations = listOfStateAfterMutations;
	}

	public List<StateAfterMutation> getListOfStateAfterMutations() {
		return listOfStateAfterMutations;
	}

	@Override
	public String toString() {
		return getCGBID() + "\t" + getID();
	}



}
