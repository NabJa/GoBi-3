package feature_extraction;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import augmentedTree.IntervalTree;
import genomicUtils.Gene;
import genomicUtils.Region;
import genomicUtils.RegionVector;
import genomicUtils.Tuple;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class SAMReader {

	public static HashMap<String, SAMRecord> lookup = new HashMap<String, SAMRecord>();
	IntervalTree<Gene> interTree = new IntervalTree<Gene>();
	Feature_writer writer;
	HashMap<Tuple<String, String>, IntervalTree<Gene>> intrTrees = new HashMap<Tuple<String, String>, IntervalTree<Gene>>();


	public SAMReader(HashMap<String, Gene> genes, String outputDestination) {
		this.interTree.addAll(genes.values());
		this.writer = new Feature_writer(outputDestination);
	}

	public void readSAM(String bamPath) {
		
		File bamF = new File(bamPath);
		SAMFileReader sam_reader = new SAMFileReader(bamF, false);
		sam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		Iterator<SAMRecord> it = sam_reader.iterator();
		
		while (it.hasNext()) {
			SAMRecord sr = it.next();

			SAMRecord other_seen = lookup.get(sr.getReadName());
			
			if (other_seen != null) {
				
				//Test from absolut read start and end if there is any invalid intergenic region:
				if (!isSuperIntergenic(other_seen, sr)) {

					String readid = sr.getReadName();
					int nsplit = getNSplit(sr, other_seen);

					if (nsplit == -1) {

						writer.writeBAMFeatureSplitIncon(readid);

					} else {

						int missmatches = getMissmatches(sr) + getMissmatches(other_seen);
						int clipping = getClipping(sr) + getClipping(other_seen);
						int gcount = 1;
						String matchGene = "Gene";
						int gdist = 1;
						boolean sense = true;
						int pcridx = 1;

						writer.writeBAMFeatures(readid, missmatches, clipping, nsplit, gcount, matchGene, gdist, sense,
								pcridx);
					}
				}
			}

			//Test if read-pair can be ignored:
			if (!can_ignore(sr)) {
				lookup.put(sr.getReadName(), sr);
			}
		}
		
		try {
			sam_reader.close();
		} catch (Exception e) {
			throw new RuntimeException("Error while closing BAM reader!", e);
		}
	}

	/**
	 * Test if there is at least one gene between read pair and no gene spanning
	 * them.
	 * 
	 * @param sr1
	 * @param sr2
	 * @return
	 */
	public boolean isSuperIntergenic(SAMRecord sr1, SAMRecord sr2) {

		
		int start = getMin(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(), sr2.getAlignmentEnd());
		int end = getMax(sr1.getAlignmentStart(), sr1.getAlignmentEnd(), sr2.getAlignmentStart(), sr2.getAlignmentEnd());
		
		int genesBetweenReads = 0;
		int readsSpanning = 0;

			ArrayList<Gene> betweenReads = interTree.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());
			for(Gene g : betweenReads) {
				if(g.getChr().equals(sr1.getReferenceName())) {
					genesBetweenReads++;					
				}
			}
			// Is there a gene between the reads?
			if (genesBetweenReads > 0) {
				
				ArrayList<Gene> genesSpanning = interTree.getIntervalsSpanning(start, end, new ArrayList<Gene>());
				for(Gene g : genesSpanning) {
					if(g.getChr().equals(sr1.getReferenceName())) {
						readsSpanning++;					
					}
				}
				// Are the reads inside of a gene?
				if (readsSpanning == 0) {
					return true;
				}
		}
		return false;
	}

	public boolean can_ignore(SAMRecord sr) {
		
		if (sr.getMateUnmappedFlag() || sr.getReadUnmappedFlag() || // is one of the reads unmapped?
				!sr.getReferenceName().equals(sr.getMateReferenceName()) || // are reads on same Chromosome?
				(sr.getReadNegativeStrandFlag() == sr.getMateNegativeStrandFlag()) || // are reads on same Strand?
				sr.getMateAlignmentStart() == 0 ||
				sr.getAlignmentStart() == 0 ||
				sr.getNotPrimaryAlignmentFlag() ||
				sr.getFlags() >= 2048 ) {
			return true;
		}

		int start1 = Math.min(sr.getAlignmentStart(), sr.getMateAlignmentStart());
		int end1 = Math.max(sr.getAlignmentStart(), sr.getMateAlignmentStart());
		int start = Math.min(start1, sr.getAlignmentEnd());
		int end = Math.max(end1, sr.getAlignmentEnd());

		
		
			ArrayList<Gene> betweenReads = interTree.getIntervalsSpannedBy(start, end, new ArrayList<Gene>());
			// Is there a gene between the reads?
			if (betweenReads.size() > 0) {
				ArrayList<Gene> genesSpanning = interTree.getIntervalsSpanning(start, end, new ArrayList<Gene>());
				// Are the reads inside of a gene?
				if (genesSpanning.size() == 0) {
					return true;
				}
		}
		return false;
	}

	public int getMissmatches(SAMRecord sr) {
		Integer nm = (Integer) sr.getAttribute("NM");
		nm = (nm != null) ? nm : (Integer) sr.getAttribute("nM");
		nm = (nm != null) ? nm : (Integer) sr.getAttribute("XM");
		return nm;
	}

	public int getClipping(SAMRecord sr) {
		int startClips = sr.getAlignmentStart() - sr.getUnclippedStart();
		int endClips = sr.getUnclippedEnd() - sr.getAlignmentEnd();
		return startClips + endClips;
	}

	public int getNSplit(SAMRecord sr, SAMRecord other_seen) {
		int nsplit = 0;
		boolean splitIntron = false;

		RegionVector readRegions = new RegionVector();
		RegionVector readRegionsRw = new RegionVector();

		for (AlignmentBlock ab : sr.getAlignmentBlocks()) {
			int ref_s = ab.getReferenceStart();
			int ref_end = ref_s + ab.getLength();

			Region readRegion = new Region(ref_s, ref_end);
			readRegions.addRegion(readRegion);
		}

		for (AlignmentBlock ab : other_seen.getAlignmentBlocks()) {

			int ref_s = ab.getReferenceStart();
			int ref_end = ref_s + ab.getLength();

			Region readRegion = new Region(ref_s, ref_end);
			readRegionsRw.addRegion(readRegion);
		}

		RegionVector readIntrons = readRegions.inverse();
		RegionVector readRwIntrons = readRegionsRw.inverse();

		for (Region intron : readIntrons.regions) {
			for (Region region : readRegionsRw.regions) {
				if (intron.getX1() == region.getX1() || intron.getX2() == region.getX2()) {
					splitIntron = true;
					break;
				}
			}
		}

		if (splitIntron) {
			return -1;
		} else {

			RegionVector uniqueIntrons = new RegionVector();
			return uniqueIntrons.regions.size();
		}
	}

	public int getMin(int v, int x, int y, int z) {
		int minVX = Math.min(v, x);
		int minYZ = Math.min(y, z);
		int gloabalMin = Math.min(minVX, minYZ); 
		return gloabalMin;
	}
	
	public int getMax(int v, int x, int y, int z) {
		int maxVX = Math.max(v, x);
		int maxYZ = Math.max(y, z);
		int globalMax = Math.max(maxVX, maxYZ); 
		return globalMax;
	}
	
	public void testIgnore(SAMRecord sr) {
		if(sr.getReadName().equals("487899")) {
			System.out.println();
			System.out.println("Stop");
			System.out.println("Stand, MateStrand " + sr.getReadNegativeStrandFlag() + " " + sr.getMateNegativeStrandFlag());
			System.out.println(sr.getReferenceName() + " " + sr.getMateReferenceName());
			System.out.println(sr.getAlignmentStart() +" " + sr.getAlignmentEnd());
			for(AlignmentBlock ab : sr.getAlignmentBlocks()) {
				System.out.println(ab.getReadStart() + " " + (ab.getReadStart() + ab.getLength()) );
			}
			

			if(sr.getFirstOfPairFlag()) {
				System.out.println("Between Reads inputs first: " + sr.getMateAlignmentStart() + " "+ sr.getAlignmentStart());
				ArrayList<Gene> betweenReads = interTree.getIntervalsSpannedBy(sr.getMateAlignmentStart(), sr.getAlignmentStart(), new ArrayList<Gene>());
				System.out.println("Between Size: " + betweenReads.size() );
				// Is there a gene between the reads?
				if (betweenReads.size() > 0) {
					ArrayList<Gene> genesSpanning = interTree.getIntervalsSpanning(sr.getMateAlignmentStart(),
							sr.getAlignmentStart(), new ArrayList<Gene>());
					// Are the reads inside of a gene?
					if (genesSpanning.size() == 0) {
						System.out.println("IGNORE");
					}
				}
			} else {
				System.out.println("Between Reads inputs second: " + sr.getAlignmentStart() + " "+ sr.getMateAlignmentStart());
				ArrayList<Gene> betweenReads = interTree.getIntervalsSpannedBy(sr.getAlignmentStart(), sr.getMateAlignmentStart(), new ArrayList<Gene>());
				System.out.println("Between Size: " + betweenReads.size() );
				// Is there a gene between the reads?
				if (betweenReads.size() > 0) {
					ArrayList<Gene> genesSpanning = interTree.getIntervalsSpanning(sr.getMateAlignmentStart(),
							sr.getAlignmentStart(), new ArrayList<Gene>());
					// Are the reads inside of a gene?
					if (genesSpanning.size() == 0) {
						System.out.println("IGNORE");
					}
				}
			}
		}
		
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
		// ArrayList<Gene> testList = testTree.getIntervalsIntersecting(100, 190, new
		// ArrayList<Gene>());
		ArrayList<Gene> testList = testTree.getIntervalsSpanning(120, 110, new ArrayList<Gene>());

		for (Gene g : testList) {
			System.out.println(g);
		}
	}

}
