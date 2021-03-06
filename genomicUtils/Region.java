package genomicUtils;

public class Region{

	public int x1;
	public int x2;
	public String regionID;
	public String geneID;
	public String bioType;

	public Region() {
		
	}
	
	public Region(int x1, int x2) {
		this.x1 = x1;
		this.x2 = x2;
	}

	 public Region(int x1, int x2, String id) {
	 this.x1 = x1;
	 this.x2 = x2;
	 this.regionID = id;
	 }
	 
	 public Region(int x1, int x2, String id, String bioType) {
	 this.x1 = x1;
	 this.x2 = x2;
	 this.regionID = id;
	 this.bioType = bioType;
	 }
	 
	 public void setRegions(int x1, int x2) {
		 this.x1 = x1;
		 this.x2 = x2;
	 }

	public String getID() {
		return regionID;
	}

	public int getX1() {
		return x1;
	}

	public int getX2() {
		return x2;
	}

	public int getLength() {
		return x2 - x1;
	}

	@Override
	public boolean equals(Object r) {
		boolean retVal = false;

        if (r instanceof Region){
        	if(((Region) r).getX1() == this.getX1() && ((Region) r).getX2() == this.getX2() ) {
        		retVal = true;
        	}
        }

     return retVal;
	}
	
	public int hashCode() {
		return ((x1 * 104723) % 104729) + ((x2 * 104717) % 104711);
	}
	
	public String toString() {
		return "" + regionID + ":" + x1 + ":" + x2;
	}
}
