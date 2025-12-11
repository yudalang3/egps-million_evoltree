package egps.module.ncov.analysis.egps.module.million.tree.treebuilder.maxparsimony;

import java.util.List;
import java.util.Objects;

import module.evoltre.mutation.IMutation4Rec;
import module.evolview.model.tree.GraphicsNode;

public class BasicOTUInfor {
	
	
	String seqAccName;
	
	
	List<IMutation4Rec> listOfMutations;
	
	GraphicsNode node = new GraphicsNode();

	int ID;
	
	private static int nextID = 0;
	
	public BasicOTUInfor() {
		this.ID = nextID;
		nextID ++;
	}

	/**
	 * @return the {@link #seqAccName}
	 */
	public String getSeqAccName() {
		return seqAccName;
	}


	/**
	 * @param seqAccName the {@link #seqAccName} to set
	 */
	public void setSeqAccName(String seqAccName) {
		this.seqAccName = seqAccName;
	}


	/**
	 * @return the {@link #listOfMutations}
	 */
	public List<IMutation4Rec> getListOfMutations() {
		return listOfMutations;
	}


	/**
	 * @param listOfMutations the {@link #listOfMutations} to set
	 */
	public void setListOfMutations(List<IMutation4Rec> listOfMutations) {
		this.listOfMutations = listOfMutations;
	}

	/**   
	 * <p>Title: hashCode</p>   
	 * <p>Description: 
	 * </p>   
	 * @return   
	 * @see java.lang.Object#hashCode()   
	 *
	 */ 
	@Override
	public int hashCode() {
		return Objects.hash(ID);
	}

	/**   
	 * <p>Title: equals</p>   
	 * <p>Description: 
	 * </p>   
	 * @param obj
	 * @return   
	 * @see java.lang.Object#equals(java.lang.Object)   
	 *
	 */ 
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (!(obj instanceof BasicOTUInfor)) {
			return false;
		}
		BasicOTUInfor other = (BasicOTUInfor) obj;
		return ID == other.ID;
	}
	
	
	public GraphicsNode getNode() {
		return node;
	}

	public int getID() {
		return ID;
	}
}
