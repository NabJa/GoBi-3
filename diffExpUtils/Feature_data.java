package diffExpUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

public class Feature_data {

	String outDestination;

	LinkedHashSet<String> allGenes = new LinkedHashSet<String>();
	LinkedHashMap<String, ArrayList<Integer>> fData = new LinkedHashMap<String, ArrayList<Integer>>();

	public Feature_data(String outputDest) {
		this.outDestination = outputDest;
	}

	public void addAllGenes(String path) {

		File geneCounts = new File(path);

		int lineCount = 0;

		try {
			BufferedReader bReader = new BufferedReader(new FileReader(geneCounts));
			String line = "";

			while ((line = bReader.readLine()) != null) {

				if (lineCount == 0)
					line = bReader.readLine();
				lineCount++;

				String geneID = line.substring(0, line.indexOf("\t"));
				allGenes.add(geneID);
			}
			bReader.close();
		} catch (Exception e) {

			e.printStackTrace();
		}

		for (String id : allGenes) {
			ArrayList<Integer> ncount = new ArrayList<Integer>();
			fData.put(id, ncount);
		}
	}

	public void addAllCounts(String path) {

		File geneCounts = new File(path);
		HashMap<String, Integer> thisData = new HashMap<String, Integer>();

		int lineCount = 0;

		try {
			BufferedReader bReader = new BufferedReader(new FileReader(geneCounts));
			String line = "";
			while ((line = bReader.readLine()) != null) {

				if (lineCount == 0)
					line = bReader.readLine();
				lineCount++;
				String[] sline = line.split("\t");
				String geneID = sline[0];
				int nreads = Integer.parseInt(sline[8]);
				thisData.put(geneID, nreads);
			}

			bReader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		for (String id : allGenes) {
			Integer thisCount = thisData.get(id);
			ArrayList<Integer> countLine = fData.get(id);

			if (thisCount == null) {
				countLine.add(0);
			} else {
				countLine.add(thisCount);
			}
		}
	}

	public void writeExprDat() {

		File exprsDat = new File(outDestination, "exprs.txt");
		BufferedWriter writer;
		FileWriter file;

		try {
			file = new FileWriter(exprsDat);
			writer = new BufferedWriter(file);

			for (String id : allGenes) {
				ArrayList<Integer> countLine = fData.get(id);
				for (int i = 0; i < countLine.size(); i++) {
					if (i < countLine.size()-1) {
						writer.write(countLine.get(i) + "\t");
					} else {
						writer.write(countLine.get(i) + "\n");
					}
				}
			}
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void writeFDat() {
		File exprsDat = new File(outDestination, "f_data.txt");
		BufferedWriter writer;
		FileWriter file;

		try {
			file = new FileWriter(exprsDat);
			writer = new BufferedWriter(file);

			for (String id : allGenes) {
				writer.write(id + "\n");
			}
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}

}
