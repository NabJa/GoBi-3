package feature_extraction;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import genomicUtils.*;

public class GTFReader {

	HashMap<String, Gene> genes = new HashMap<String, Gene>();

	public HashMap<String, Gene> getGenes(){
		return genes;
	}
	
	public void readExon(String path) {

		BufferedReader reader = null;

		try {
			File file = new File(path);
			reader = new BufferedReader(new FileReader(file));

			String rline = "";
			
			while ((rline = reader.readLine()) != null) {

				if (rline.indexOf('#') != 0) {

					int firstTab = rline.indexOf('\t');
					int secondTab = rline.indexOf('\t', firstTab + 1);
					
					int thirdTab = rline.indexOf('\t', secondTab + 1);
					int fourthTab = rline.indexOf('\t', thirdTab + 1);
					int fifthTab = rline.indexOf('\t', fourthTab + 1);
					int sixthTab = rline.indexOf('\t', fifthTab + 1);
					int seventhTab = rline.indexOf('\t', sixthTab + 1);

					if (rline.substring(secondTab + 1, secondTab + 5).toLowerCase().equals("exon")) {
											
						String chr = rline.substring(0, firstTab);
						String bioType = rline.substring(firstTab + 1, secondTab);
						int start = Integer.parseInt(rline.substring(thirdTab + 1, fourthTab));
						int end = Integer.parseInt(rline.substring(fourthTab + 1, fifthTab));
						String strand = rline.substring(sixthTab + 1, seventhTab);
						
						String[] line = rline.split(";");
						
						String transID = transcriptIDSearch(line);
						String geneID = geneIDSearch(line);
//						String geneName = geneNameSearch(line);
						String proteinID = proteinIDSearch(line);

						Gene gene = new Gene();
						gene.setGene(geneID, chr, strand, start, end, bioType);
						Gene newGene = genes.putIfAbsent(geneID, gene);

						Region cds = new Region(start, end, proteinID, bioType);
						RegionVector transcript = new RegionVector(transID);
						
						
//						System.out.println(chr + "\t" + bioType + "\t" + start  + "\t" + end  + "\t" + strand  + "\t" + transID  + "\t" + geneID + "\t" + proteinID);
						
						if (newGene == null) // means this gene is new
						{
							transcript.addRegion(cds);
							gene.transcripts.put(transID, transcript);

						} else // gene already exists
						{
							Gene correspondingGene = genes.get(geneID);
							RegionVector newTrans = correspondingGene.transcripts.putIfAbsent(transID, transcript);

							if (newTrans == null) // means this Trans is new
							{
								transcript.addRegion(cds);
								correspondingGene.transcripts.put(transID, transcript);
							} else // Trans already exists
							{
								genes.get(geneID).transcripts.get(transID).addRegion(cds);
							}
						}
					}
				}
			}

		} catch (Exception e) {
			throw new RuntimeException("got error while reading gtf.", e);
		} finally {

			try {
				reader.close();
			} catch (Exception e) {
				throw new RuntimeException("got error while closing gtf.", e);
			}

		}
	}

	public String transcriptIDSearch(String[] line) {
		String targetID = "";
		for (int i = 0; i < line.length; i++) {
			if (line[i].indexOf('t') == 1 && line[i].indexOf('d') == 13) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public String proteinIDSearch(String[] line) {
		String targetID = "";
		for (int i = 0; i < line.length; i++) {
			if (line[i].indexOf('p') == 1 && line[i].indexOf('d') == 10) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public String geneNameSearch(String[] line) {
		String targetID = "";
		for (int i = 0; i < line.length; i++) {
			if (line[i].indexOf('g') == 1 && line[i].indexOf('m') == 8) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public String geneIDSearch(String[] line) {
		String targetID = "";
		for (int i = 0; i < line.length; i++) {
			line[i] = line[i].trim();
			if (line[i].indexOf('g') == 0 && line[i].indexOf('d') == 6) {
				targetID = line[i].substring(line[i].indexOf('\"') + 1, line[i].length() - 1);
				break;
			}
		}
		return targetID;
	}

	public boolean geneAnnotated(String geneID) {
		boolean bol = false;
		for (String ids : genes.keySet()) {
			if (geneID == ids) {
				bol = true;
				break;
			}
		}
		return bol;
	}


}

