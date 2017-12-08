package genomicUtils;

public class TupleComperator implements Comperator<Tuple<String, String>> {

	@Override
	public int compare(Tuple<String, String> o1, Tuple<String, String> o2) {

		int diff = 1;
		if (o1.getFirst().equals(o2.getFirst()) && o1.getSecond().equals(o2.getSecond())) {
			diff = 0;
		}
		return diff;
	}
}
