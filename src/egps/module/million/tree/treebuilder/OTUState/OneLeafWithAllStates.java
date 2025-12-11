package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.OTUState;

import module.parsimonytre.algo.StateAfterMutation;

import java.util.List;

public class OneLeafWithAllStates {


    List<StateAfterMutation> allStates;


    public void setAllStates(List<StateAfterMutation> allStates) {
        this.allStates = allStates;
    }

    public List<StateAfterMutation> getAllStates() {

        return allStates;
    }
}
