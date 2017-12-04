package feature_extraction;

import java.util.ArrayList;

import genomicUtils.*;



public class Feature_extraction {
	
		ArrayList<Triplet<String, String, String>> bams = new ArrayList<Triplet<String, String, String>>();

		
		public static void main(String args[]) {
		
//			String gtf = "/home/j/jabareen/Desktop/GoBi/GoBi-3/BamFeatures/complete_bams/Homo_sapiens.GRCh37.75.gtf";
//			String bam = "/home/j/jabareen/Desktop/GoBi/GoBi-3/BamFeatures/complete_bams/hes_star.bam";

			String gtf = "/home/j/jabareen/Desktop/GoBi/GoBi-3/BamFeatures/complete_bams/Saccharomyces_cerevisiae.R64-1-1.75.gtf";
			String bam = "/home/j/jabareen/Desktop/GoBi/GoBi-3/BamFeatures/complete_bams/nookaew_cm.bam";

			
			GTFReader gtfReader = new GTFReader();
			System.out.println("Start reading GTF");
			gtfReader.readExon(gtf);
						
//			for(Gene g : gtfReader.genes.values()) {
//				System.out.println("Gene: " + g.geneID + " Start: " + g.start + " End: " + g.end + " Length: " + (g.end-g.start));
//			}
			
			SAMReader samReader = new SAMReader(gtfReader.genes);
			System.out.println("Start reading SAM");
			samReader.readSAM(bam);
		
		
		}
}
 	