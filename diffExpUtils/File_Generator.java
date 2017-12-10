package diffExpUtils;

import java.io.File;

public class File_Generator {

	static String countFiles;
	static String labels;
	static String output;

	public static void main(String args[]) {

		readArguments(args);

		File countF = new File(countFiles);

		Countfiles_Reader countReader = new Countfiles_Reader(countF);
		countReader.readCountFiles(output);

	}

	public static void readArguments(String args[]) {

		for (int i = 0; i < args.length; i++) {

			switch (args[i]) {
			case "-countfiles":
				countFiles = args[i + 1];
				System.out.println("Countfiles: " + countFiles);
				break;

			case "-labels":
				labels = args[i + 1];
				System.out.println("Labels: " + labels);
				break;

			case "-o":
				output = args[i + 1];
				System.out.println("Output: " + output);
				break;

			default:
				break;
			}
		}
		System.out.println();
	}
}
