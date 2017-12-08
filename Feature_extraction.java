package feature_extraction;

import java.util.ArrayList;

import genomicUtils.*;

public class Feature_extraction {

	ArrayList<Triplet<String, String, String>> bams = new ArrayList<Triplet<String, String, String>>();

	static String gtf;
	static String bam;
	static String outDest;
	static boolean frstrand;
	
	public static void main(String args[]) {

		readArguments(args);
		
		GTFReader gtfReader = new GTFReader();
		System.out.println("Start reading GTF");
		gtfReader.readExon(gtf);

		SAMReader samReader = new SAMReader(gtfReader.genes, outDest);
		System.out.println("Start reading SAM");
		samReader.readSAM(bam);
		samReader.writer.closeBAMFeatures();

	}

	public static void readArguments(String args[]) {

		String noInp = "Error while reading input instructions!!";

		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-gtf":
				gtf = args[i + 1];
				i++;
				System.out.println("GTF: " + gtf);
				break;

			case "-bam":
				bam = args[i + 1];
				i++;
				System.out.println("BAM: " + bam);
				break;
				
			case "-o":
				outDest = args[i + 1];
				i++;
				System.out.println("Output Destination: " + outDest);
				break;
				
			case "-frstrand":
				frstrand = Boolean.parseBoolean(args[i + 1]);
				i++;
				System.out.println("Fragment strand: " + frstrand);
				break;
			
			default:
				System.out.println(noInp);
			}
		}

	}
}
