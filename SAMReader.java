package feature_extraction;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import augmentedTree.IntervalTree;
import genomicUtils.Gene;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class SAMReader {

	public static HashMap<String, SAMRecord> lookup = new HashMap<String, SAMRecord>();
	IntervalTree<Gene> interTree = new IntervalTree<Gene>();

	public SAMReader(HashMap<String, Gene> genes) {
		this.interTree.addAll(genes.values());
	}

	public void readSAM(String bamPath) {

		File bamF = new File(bamPath);
		SAMFileReader sam_reader = new SAMFileReader(bamF, false);
		sam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		Iterator<SAMRecord> it = sam_reader.iterator();

	
		while (it.hasNext()) {
			SAMRecord sr = it.next();

			if (can_ignore(sr)) {
				lookup.put(sr.getReadName(), sr);
			}

			SAMRecord other_seen = lookup.get(sr.getReadName());
			if (other_seen != null) {

			}
		}

		try {
			sam_reader.close();
		} catch (Exception e) {
			throw new RuntimeException("Error while closing BAM reader!", e);
		}

	}

	public boolean can_ignore(SAMRecord sr) {
		boolean re = false;

		ArrayList<Gene> intersectingGenes = interTree.getIntervalsIntersecting(sr.getAlignmentStart(), sr.getMateAlignmentStart(), new ArrayList<Gene>());

		if (sr.getMateUnmappedFlag() || sr.getMateReferenceName() == null || sr.getMateAlignmentStart() == 0 || intersectingGenes.size() < 1) {
			re = true;

		}

		return re;
	}

	public void testIntervallTree() {

		HashMap<String, Gene> testHash = new HashMap<String, Gene>();

		Gene gene1 = new Gene(10, 20);
		Gene gene2 = new Gene(30, 40);
		Gene gene3 = new Gene(100, 150);
		Gene gene4 = new Gene(180, 220);

		testHash.put("gene1", gene1);
		testHash.put("gene2", gene2);
		testHash.put("gene3", gene3);
		testHash.put("gene4", gene4);

		IntervalTree<Gene> testTree = new IntervalTree<Gene>();
		testTree.addAll(testHash.values());
		ArrayList<Gene> testList = testTree.getIntervalsIntersecting(40, 100, new ArrayList<Gene>());

		for (Gene g : testList) {
			System.out.println(g);
		}
	}

}
