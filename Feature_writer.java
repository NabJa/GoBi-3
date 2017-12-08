package feature_extraction;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class Feature_writer {

	String outputDestination;
	FileWriter file;
	BufferedWriter writer;

	public Feature_writer(String outputDestination) {
		this.outputDestination = outputDestination;

		try {
			this.file = new FileWriter(outputDestination);
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while creating BufferedWriter for Feature_writer.", e);
		}
	}

	public void writeBAMFeatures(String readid, int mm, int clipping, int nsplit, int gcount, String matchGene,
			int gdist, boolean sense, int pcridx) {
		try {

			writer.write(readid + "\t" + "mm:" + mm + "\t" + "clipping:" + clipping + "\t" + "nsplit:" + nsplit + "\t"
					+ "gcount:" + gcount + "\t" + matchGene + "\t" + "gdist:" + gdist + "\t" + sense + "\t" + "pcrindex"
					+ pcridx);
			writer.newLine();
		} catch (Exception e) {
			throw new RuntimeException("Got error while writing Feature writer." + e);
		}
	}

	public void writeBAMFeatureSplitIncon(String readid) {
		try {

			writer.write(readid + "\t" + "split-inconsistent:true");
			writer.newLine();
			
		} catch (Exception e) {
			throw new RuntimeException("Got error while writing Feature writer." + e);
		}
	}

	public void closeBAMFeatures() {
		try {
			writer.close();
		} catch (Exception e) {
			throw new RuntimeException("Got error while closing Feature writer", e);
		}
	}

}
