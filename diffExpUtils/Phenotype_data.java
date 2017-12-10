package diffExpUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class Phenotype_data {

	String outDestination;
	
	ArrayList<String> phenotype = new ArrayList<String>();
	ArrayList<Integer> feature = new ArrayList<Integer>();
	HashMap<String, Integer> seen_feature = new HashMap<String, Integer>();
	int feat = -1;

	public Phenotype_data(String outDestination) {
	this.outDestination = outDestination;
	}

	public void addPhenotype(String condition, String replicate) {
		String phen = condition + "." + replicate;
		phenotype.add(phen);
	}

	public void addFeature(String new_feat) {

		if (seen_feature.containsKey(new_feat)) {
			int seen_feat = seen_feature.get(new_feat);
			feature.add(seen_feat);
		} else {
			feat++;
			seen_feature.put(new_feat, feat);
			feature.add(feat);
		}
	}

	public void writePDat() {
		File exprsDat = new File(outDestination, "p_data.txt");
		BufferedWriter writer;
		FileWriter file;

		try {
			file = new FileWriter(exprsDat);
			writer = new BufferedWriter(file);

			for (int i = 0; i < phenotype.size(); i++) {
			writer.write(phenotype.get(i) + "\t" + feature.get(i) + "\n");
			}
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	public void printPDat() {
		for (int i = 0; i < phenotype.size(); i++) {
			System.out.print(phenotype.get(i) + "\t" +"\t" + feature.get(i));
			System.out.println();
		}
	}

}
