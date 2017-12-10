package diffExpUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class Countfiles_Reader {
	
	File countFile;
	
	public Countfiles_Reader(File countFile) {
		this.countFile = countFile;
	}
	
	public void readCountFiles(String output) {
		
		Phenotype_data pData = new Phenotype_data(output);
		Feature_data fData = new Feature_data(output);
		
		try {
			
			BufferedReader bReader = new BufferedReader(new FileReader(countFile));
			
			String line = "";
			
			while((line = bReader.readLine()) != null) {
				
				String[] sline = line.split("\t");
				
				pData.addFeature(sline[0]);
				pData.addPhenotype(sline[0], sline[2]);

//				System.out.println(sline[1]);
				
				fData.addAllGenes(sline[1]);
				
			}
			bReader.close();
			
			BufferedReader cReader = new BufferedReader(new FileReader(countFile));
			String cline = "";

			while((cline =cReader.readLine()) != null) {
				
				String[] sline = cline.split("\t");
				
				fData.addAllCounts(sline[1]);
				
			}
			cReader.close();
			
			pData.writePDat();
			fData.writeFDat();
			fData.writeExprDat();
					
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
