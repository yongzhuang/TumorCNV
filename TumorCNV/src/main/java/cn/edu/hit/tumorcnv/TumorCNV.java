package cn.edu.hit.tumorcnv;

import cn.edu.hit.tumorcnv.Preprocessing;
import java.io.File;
import java.io.IOException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;

/**
 * @author Yongzhuang Liu
 *
 */
public class TumorCNV {

	private static Logger logger = Logger.getLogger(TumorCNV.class);

	public static void main(String[] args) throws IOException {
		String usage = "\nTumorCNV-0.1.0\n";
		usage = usage + "\nUsage: java -jar TumorCNV.jar <COMMAND> [OPTIONS]\n\n";
		usage = usage + "COMMANDS:\n" + "\tpreprocess\textract information\n"
				+ "\tcall\t\tcall germline and somatic CNVs \n";
		if (args.length > 0) {
			if (args[0].equals("preprocess") || args[0].equals("call")) {
				run(args[0], args);
			} else {
				logger.error("Command is not recognized!\n" + usage);
			}
		} else {
			System.out.println(usage);
		}
	}

	private static void run(String cmd, String[] args) throws IOException {
		long start = System.currentTimeMillis();
		CommandLineParser parser = new PosixParser();
		CommandLine commandLine = null;
		Options options = TumorCNV.createOptions(cmd);
		try {
			if (options != null) {
				commandLine = parser.parse(options, args);
			}
		} catch (ParseException parseException) {
			logger.error("Invalid command line parameters!");
		}
		if (cmd.equals("preprocess")) {
			if (isValidated(commandLine, "preprocess")) {
				String normalVCFFile = commandLine.getOptionValue("normalVCFFile");
				String normalBAMFile = commandLine.getOptionValue("normalBAMFile");
				String tumorBAMFile = commandLine.getOptionValue("tumorBAMFile");
				String referenceSequenceFile = commandLine.getOptionValue("referenceSequenceFile");
				String mappabilityFile = commandLine.getOptionValue("mappabilityFile");
				String outputPrefix = commandLine.getOptionValue("outputPrefix");
				int minMappingQuality = Integer.parseInt(commandLine.getOptionValue("minMappingQuality", "1"));
				int minBaseQuality = Integer.parseInt(commandLine.getOptionValue("minBaseQuality", "20"));
				int windowSize = Integer.parseInt(commandLine.getOptionValue("windowSize", "500"));
				Preprocessing preprocessing = new Preprocessing(normalVCFFile, normalBAMFile, tumorBAMFile,
						referenceSequenceFile, mappabilityFile, outputPrefix, minMappingQuality, minBaseQuality,
						windowSize);
				preprocessing.run();
			} else {
				printHelp("preprocess");
				return;
			}
		}
		if (cmd.equals("call")) {
			if (isValidated(commandLine, "call")) {
				String inputBAFFile = commandLine.getOptionValue("afFile");
				String inputRDFile = commandLine.getOptionValue("rdFile");
				String mappabilityFile = commandLine.getOptionValue("mappabilityFile");
				String outputPrefix = commandLine.getOptionValue("outputPrefix");
				String intervalFile = commandLine.getOptionValue("exclude");
				double minMappability = Double.parseDouble(commandLine.getOptionValue("minMappability", "0.30"));
				int minDistance = Integer.parseInt(commandLine.getOptionValue("minDistance", "10000"));
				double transitionProb = Double.parseDouble(commandLine.getOptionValue("transitionProb", "0.00001"));
				double outlier = Double.parseDouble(commandLine.getOptionValue("outlier", "0.1"));
				double purity = Double.parseDouble(commandLine.getOptionValue("purity", "1.0"));
				int ploidy = Integer.parseInt(commandLine.getOptionValue("ploidy", "2"));
				int numOfThreads = Integer.parseInt(commandLine.getOptionValue("nt", "1"));
				new Calling(inputBAFFile, inputRDFile, mappabilityFile, outputPrefix, intervalFile, transitionProb,
						outlier, purity, ploidy, minMappability, minDistance).runMultiThreading(numOfThreads);
			} else {
				printHelp("call");
				return;
			}
		}
		long end = System.currentTimeMillis();
		logger.info("Total running time is " + (end - start) / 1000 + " seconds");
		logger.info("Done!");
	}

	private static Options createOptions(String cmd) {
		Options options = new Options();
		if (cmd.equals("preprocess")) {
			options.addOption(OptionBuilder.withLongOpt("normalVCFFile").withDescription("normal vcf file (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("normalBAMFile").withDescription("normal bam file (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("tumorBAMFile").withDescription("tumor bam file (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("referenceSequenceFile")
					.withDescription("reference genome file (required)").hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("mappabilityFile")
					.withDescription("mappability file (required)").hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("outputPrefix")
					.withDescription("perfix of output file (required)").hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("minMappingQuality")
					.withDescription("minumum mapping quality (optional, default 1)").hasArg().withArgName("INT")
					.create(""));
			options.addOption(OptionBuilder.withLongOpt("minBaseQuality")
					.withDescription("minumum base quality (optional, default 20)").hasArg().withArgName("INT")
					.create(""));
			options.addOption(OptionBuilder.withLongOpt("windowSize")
					.withDescription("window size (optional, default 500)").hasArg().withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("purity").withDescription("tumor purity").hasArg()
					.withArgName("INT").create());
			return options;
		} else if (cmd.equals("call")) {
			options.addOption(OptionBuilder.withLongOpt("afFile").withDescription("allele frequency file (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("rdFile").withDescription("the RD file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("outputPrefix").withDescription("the output prefix (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("mappabilityFile")
					.withDescription("mappability file (required)").hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("exclude").withDescription("exclude regions").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("transitionProb").withDescription("transitionProb").hasArg()
					.withArgName("DOUBLE").create());
			options.addOption(OptionBuilder.withLongOpt("minMappability").withDescription("minMappability").hasArg()
					.withArgName("DOUBLE").create());
			options.addOption(OptionBuilder.withLongOpt("minDisatance").withDescription("minDisatance").hasArg()
					.withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("outlier").withDescription("outlier").hasArg()
					.withArgName("DOUBLE").create());
			options.addOption(OptionBuilder.withLongOpt("purity").withDescription("purity").hasArg()
					.withArgName("DOUBLE").create());
			options.addOption(
					OptionBuilder.withLongOpt("ploidy").withDescription("ploidy").hasArg().withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("nt").withDescription("number of threads (optional, default 1)")
					.hasArg().withArgName("INT").create());
			return options;
		} else {
			return null;
		}
	}

	private static boolean isValidated(CommandLine line, String cmd) {
		boolean tag = true;
		if (cmd.equals("preprocess")) {
			if (!line.hasOption("normalBAMFile") || !(new File(line.getOptionValue("normalBAMFile")).isFile())) {
				logger.error("The normal sample's BAM file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("tumorBAMFile") || !(new File(line.getOptionValue("normalBAMFile")).isFile())) {
				logger.error("The tumor sample's BAM file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("referenceSequenceFile")
					|| !(new File(line.getOptionValue("referenceSequenceFile")).isFile())) {
				logger.error("The reference genome file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("mappabilityFile") || !(new File(line.getOptionValue("mappabilityFile")).isFile())) {
				logger.error("The mappability file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("outputPrefix")) {
				logger.error("The output prefix is not correctly specified!");
				tag = false;
			}
		}
		if (cmd.equals("call")) {
			if (!line.hasOption("rdFile") || !(new File(line.getOptionValue("rdFile")).isFile())) {
				logger.error("The read depth file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("mappabilityFile") || !(new File(line.getOptionValue("mappabilityFile")).isFile())) {
				logger.error("The mappability file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("outputPrefix")) {
				logger.error("The output prefix is not correctly specified!");
				tag = false;
			}
		}
		return tag;
	}

	private static void printHelp(String command) {
		System.out.println();
		String usage1 = "usage: java -jar TumorCNV.jar " + command + " [OPTIONS]\n\n"
				+ "-referenceSequenceFile\t<FILE>\treference genome file (required)\n"
				+ "-normalVCFFile\t<FILE>\tnormal sample's vcf file (optional) \n"
				+ "-normalBAMFile\t<FILE>\tnormal sample's bam file (required)\n"
				+ "-tumorBAMFile\t<FILE>\ttumor sample's bam file (required)\n"
				+ "-mappabilityFile\t<FILE>\tmappability file (required)\n"
				+ "-outputPrefix\t<FILE>\t prefix of output file (required)\n"
				+ "-windowSize\t<INT>\twindow size (optional, default 500)\n"
				+ "-minMappingQuality\t<INT>\tminimum mapping quality (optional, default 1)\n"
				+ "-minBaseQuality\t<INT>\tminimum base quality (optional, default 20)\n";
		String usage2 = "usage: java -jar TumorCNV.jar " + command + " [OPTIONS]\n\n"
				+ "-rdFile\t<FILE>\tread depth file (required)\n"
				+ "-afFile\t<FILE>\tallele frequency file (optional)\n"
				+ "-mappabilityFile\t<FILE>\tmappability file (required)\n"
				+ "-outputPrefix\t<FILE>\tprefix of toutput file (required)\n" + "-exclude\t<FILE>\texclude regions\n"
				+ "-transitionProb\t<FLOAT>\ttransition probability of different states (optional, default 0.00001)\n"
				+ "-minMappability\t<FLOAT>\tminimum mappability of window (optional, default 0.3)\n"
				+ "-minDisatance\t<INT>\tminimum distance to merge adjacent CNVs (optional, default 10000)\n"
				+ "-purity\t<FLOAT>\ttumor purity (optional, default 1.0)\n"
				+ "-ploidy\t<INT>\ttumor ploidy (optional, default 2)\n"
				+ "-outlier\t<FLOAT>\tthe percentage of outliers (optional, default 0.1)\n"
				+ "-nt\t<INT>\tnumber of threads (optional, default 1)\n";
		if (command.equals("preprocess")) {
			System.out.println(usage1);
		}
		if (command.equals("call")) {
			System.out.println(usage2);
		}
	}
}
