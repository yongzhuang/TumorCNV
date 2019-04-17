package cn.edu.hit.tumorcnv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import rcaller.RCaller;
import rcaller.RCode;

/**
 * @author Yongzhuang Liu
 *
 */
public class Calling {

	private static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(Calling.class);

	private String rdFile;
	private String outputFile;
	private String mappabilityFile;
	private String intervalFile;
	private String afFile;
	private double minMappability;
	private int minDistance;
	private double outlier;
	private double tau;
	private double alpha;
	private double p;

	private static final int MIN_GC_BIN_SIZE = 1000;
	private static final int MIN_GC = 30;
	private static final int MAX_GC = 70;

	public Calling(String afFile, String rdFile, String mappabilityFile, String outputFile, String intervalFile,
			double transitionProb, double outlier, double purity, int ploidy, double minMappability, int minDistance) {
		this.afFile = afFile;
		this.rdFile = rdFile;
		this.mappabilityFile = mappabilityFile;
		this.outputFile = outputFile;
		this.outlier = outlier;
		this.alpha = purity;
		this.tau = ploidy;
		this.minMappability = minMappability;
		this.minDistance = minDistance;
		this.p = transitionProb;
		this.intervalFile = intervalFile;
	}

	public void runMultiThreading(int numOfThreads) throws IOException {
		logger.info("Checking R environment and packages ......");
		try {
			Process pid = Runtime.getRuntime().exec("which Rscript");
			BufferedReader runTimeReader = new BufferedReader(new InputStreamReader(pid.getInputStream()));
			String R_HOME = runTimeReader.readLine().trim();
			if (R_HOME.equals("")) {
				logger.error("Rscript exectuable is not set in the PATH environment variable!");
				return;
			}
			if (!checkInstalledPackages(R_HOME)) {
				return;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		logger.info("Loading data ......");
		List<Observation> observationList = null;
		if (this.afFile == null) {
			observationList = getObservationList(rdFile);
		} else {
			observationList = getObservationList(rdFile, afFile);
		}
		List<List<Observation>> rawObservationList = getObservationListByChrom(observationList);
		List<List<Observation>> observationListByChrom;
		if (this.intervalFile != null) {
			IntervalReader intervalReader = new IntervalReader(intervalFile);
			List<Interval> intervalList = intervalReader.getIntervalList();
			observationListByChrom = excludeRegions(rawObservationList, intervalList);
		} else {
			observationListByChrom = rawObservationList;
		}
		File file = new File(outputFile);
		if (file.exists()) {
			try {
				file.delete();
				file.createNewFile();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		ThreadPoolExecutor threadPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(numOfThreads);
		for (int i = 0; i < observationListByChrom.size(); i++) {
			List<Observation> observations = observationListByChrom.get(i);
			String chrom = observations.get(0).getChrom();
			logger.info("Chromosome " + chrom + " ......");
	        if(chrom.equals("X")||chrom.equals("Y")) {
	        	break;
	        }
			List<Observation> observationList1 = getFilteredReadDepthList(observations, MIN_GC, MAX_GC, minMappability);
			double rho = 0;
			NBModel[] nbModels = getNBModels(observationList1);
			Emission emission = new Emission(nbModels[0], nbModels[1], alpha, tau, rho);
			List<Observation> observationList2 = null;
			if (this.afFile != null) {
				double pvalue = testOverdispersion(observationList1);
				double[] parameters = getDispersion(observationList1);
				if (pvalue < 0.05) {
					rho = Math.max(parameters[2], 0.001);
				}
				observationList2 = getFilteredAlleleFrequencyList(observationList1, parameters[0], parameters[1]);
			} else {
				observationList2 = observationList1;
			}
			if (observationList2.size() > 0) {
				double[] pi = getPI();
				int D = observationList2.get(0).getEnd() - observationList2.get(0).getStart() + 1;
				Transition transition = (new Transition(D, p));
				Viterbi viterbi2 = (new Viterbi(observationList2, pi, transition, emission));
				threadPool.execute(new SingleThreadCall(viterbi2, outputFile, minDistance, tau));
			}
		}
		threadPool.shutdown();
		try {
			threadPool.awaitTermination(30, TimeUnit.DAYS);
			boolean loop = true;
			long memory = 0;
			long maxTotalMemory = 0;
			do {
				long totalMemory = Runtime.getRuntime().totalMemory();
				long freeMemory = Runtime.getRuntime().freeMemory();
				long tmp = (totalMemory - freeMemory) / (1024 * 1024);
				if (tmp > memory) {
					memory = tmp;
				}
				if (totalMemory > maxTotalMemory) {
					maxTotalMemory = totalMemory;
				}
				loop = !threadPool.awaitTermination(1, TimeUnit.MINUTES);
			} while (loop);
			threadPool.awaitTermination(30, TimeUnit.DAYS);
		} catch (InterruptedException ex) {
			ex.printStackTrace();
		}

		List<CNVRecord> records = new ArrayList<CNVRecord>();
		try {
			BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(this.outputFile)));
			String line;
			while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
				String[] items = line.split("\t");
				String chrom = items[0];
				int start = Integer.parseInt(items[1]);
				int end = Integer.parseInt(items[2]);
				int[] states = new int[2];
				states[0] = Integer.parseInt(items[3]);
				states[1] = Integer.parseInt(items[4]);
				CNVTYPE type = parseCNVType(items[5]);
				double mappability = Double.parseDouble(items[6]);
				records.add(new CNVRecord(chrom, start, end, states, type, mappability));
			}
			bufferedReader.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		Collections.sort(records);
		file.delete();
		try {
			PrintWriter writer1 = new PrintWriter(this.outputFile + ".Germline");
			PrintWriter writer2 = new PrintWriter(this.outputFile + ".Somatic");
			writer1.write("CHROM\tSTART\tEND\tTYPE\tMAPPABILITY\tLENGTH\n");
			writer2.write("CHROM\tSTART\tEND\tTYPE\tMAPPABILITY\tLENGTH\n");
			for (CNVRecord record : records) {
				double mappability = getMappability(mappabilityFile, record.getChrom(), record.getStart(),
						record.getEnd());
				CNVTYPE type = record.getType();
				if (type.equals(CNVTYPE.GERMLINE_DEL)) {
					writer1.write(record.getChrom() + "\t" + record.getStart() + "\t" + record.getEnd() + "\t" + "DEL\t"
							+ mappability + "\t" + record.getLength() + "\n");
				}
				if (type.equals(CNVTYPE.GERMLINE_DUP)) {
					writer1.write(record.getChrom() + "\t" + record.getStart() + "\t" + record.getEnd() + "\t" + "DUP\t"
							+ mappability + "\t" + record.getLength() + "\n");
				}
				if (type.equals(CNVTYPE.SOMATIC_GAIN)) {
					writer2.write(record.getChrom() + "\t" + record.getStart() + "\t" + record.getEnd() + "\t"
							+ "GAIN\t" + mappability + "\t" + record.getLength() + "\n");
				}
				if (type.equals(CNVTYPE.SOMATIC_LOSS)) {
					writer2.write(record.getChrom() + "\t" + record.getStart() + "\t" + record.getEnd() + "\t"
							+ "LOSS\t" + mappability + "\t" + record.getLength() + "\n");
				}
			}
			writer1.close();
			writer2.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private List<Observation> getFilteredReadDepthList(List<Observation> observations, int minGC, int maxGC,
			double minMappability) {
		List<Observation> newObservations = new ArrayList();
		for (int j = 0; j < observations.size(); j++) {
			Observation tmp = observations.get(j);
			if (tmp.getGC() >= minGC && tmp.getGC() <= maxGC && tmp.getMappability() >= minMappability) {
				newObservations.add(tmp);
			}
		}
		return newObservations;
	}

	private List<Observation> getFilteredAlleleFrequencyList(List<Observation> observations, double minAF,
			double maxAF) {
		List<Observation> newObservations = new ArrayList();
		for (int j = 0; j < observations.size(); j++) {
			Observation tmp = observations.get(j);
			if (tmp.getNormalAlleleCountList() != null) {
				List<int[]> normalAlleleCountList = tmp.getNormalAlleleCountList();
				List<int[]> tumorAlleleCountList = tmp.getTumorAlleleCountList();
				List<int[]> newNormalAlleleCountList = new ArrayList();
				List<int[]> newTumorAlleleCountList = new ArrayList();
				for (int i = 0; i < normalAlleleCountList.size(); i++) {
					double alleleFrequency = (double) normalAlleleCountList.get(i)[0]
							/ (double) normalAlleleCountList.get(i)[0] + (double) normalAlleleCountList.get(i)[1];
					if (alleleFrequency > minAF && alleleFrequency < maxAF) {
						newNormalAlleleCountList.add(normalAlleleCountList.get(i));
						newTumorAlleleCountList.add(tumorAlleleCountList.get(i));
					}
				}
				tmp.setNormalAlleleCountList(newNormalAlleleCountList);
				tmp.setTumorAlleleCountList(newTumorAlleleCountList);
			}
			newObservations.add(tmp);
		}
		return newObservations;
	}

	private double getMappability(String mappabilityFile, String chrom, int start, int end) {
		if (!chrom.startsWith("chr")) {
			chrom = "chr" + chrom;
		}
		try {
			BBFileReader reader = new BBFileReader(mappabilityFile);
			BigWigIterator iterator = reader.getBigWigIterator(chrom, start - 1, chrom, end, false);
			double total = 0;
			while (iterator.hasNext()) {
				WigItem item = iterator.next();
				double value = item.getWigValue();
				int startBase = item.getStartBase();
				int endBase = item.getEndBase();
				if (startBase < start - 1 && endBase > end) {
					total += value * (end - start + 1);
				} else if (startBase < start - 1 && endBase <= end) {
					total += value * (endBase - start + 1);
				} else if (startBase >= start - 1 && endBase > end) {
					total += value * (end - startBase);
				} else {
					total += value * (endBase - startBase);
				}
			}
			reader.close();
			return total / (end - start + 1);
		} catch (IOException ex) {
			ex.printStackTrace();
			return 0;
		}
	}

	public List<Observation> getObservationList(String readDepthFile) throws IOException {
		List<Observation> observationList = new ArrayList();
		BufferedReader rdFileReader = new BufferedReader(new FileReader(readDepthFile));
		String rdLine = rdFileReader.readLine();
		while ((rdLine = rdFileReader.readLine()) != null) {
			String[] record = rdLine.split("\t");
			String chrom = record[0];
			int start = Integer.parseInt(record[1]);
			int end = Integer.parseInt(record[2]);
			int[] rd = new int[] { Integer.parseInt(record[3]), Integer.parseInt(record[4]) };
			int gc = Integer.parseInt(record[5]);
			double mappability = Double.parseDouble(record[6]);
			Observation observation = new Observation(chrom, start, end, gc, mappability, null, null, rd);
			observationList.add(observation);
		}
		rdFileReader.close();
		return observationList;
	}

	public List<Observation> getObservationList(String readDepthFile, String alleleFrequencyFile) throws IOException {
		List<Observation> observationList = new ArrayList();
		List<Observation> rawObservationList = getObservationList(readDepthFile);
		String currentChrom = null;
		List<Observation> currentChromObservationList = null;
		BufferedReader afFileReader = new BufferedReader(new FileReader(alleleFrequencyFile));
		String afLine = afFileReader.readLine();
		while ((afLine = afFileReader.readLine()) != null) {
			String[] record = afLine.split("\t");
			String chrom = record[0];
			int pos = Integer.parseInt(record[1]);
			int normalRefCount = Integer.parseInt(record[2]);
			int normalAltCount = Integer.parseInt(record[3]);
			int tumorRefCount = Integer.parseInt(record[4]);
			int tumorAltCount = Integer.parseInt(record[5]);
			if (currentChrom == null) {
				currentChrom = chrom;
				currentChromObservationList = getObservationListOfChrom(rawObservationList, chrom);
			}
			if (currentChrom != null && !currentChrom.equals(chrom)) {
				observationList.addAll(currentChromObservationList);
				currentChrom = chrom;
				currentChromObservationList = getObservationListOfChrom(rawObservationList, chrom);
				if (currentChromObservationList == null) {
					break;
				}
			}
			int windowSize = rawObservationList.get(0).getEnd() - rawObservationList.get(0).getStart() + 1;
			int index = (pos - 1) / windowSize;
			int startIndex = (currentChromObservationList.get(0).getStart() - 1) / windowSize;
			if (index - startIndex < currentChromObservationList.size()) {
				currentChromObservationList.get(index - startIndex)
						.addNormalAlleleCountList(new int[] { normalRefCount, normalAltCount });
				currentChromObservationList.get(index - startIndex)
						.addTumorAlleleCountList(new int[] { tumorRefCount, tumorAltCount });
			}
		}
		if (currentChromObservationList != null) {
			observationList.addAll(currentChromObservationList);
		}
		afFileReader.close();
		return observationList;
	}

	private List<List<Observation>> getObservationListByChrom(List<Observation> observationList) {
		List<List<Observation>> observations = new ArrayList();
		List<Observation> tmp = null;
		String lastChrom = null;
		int lastEnd = 0;
		for (Observation observation : observationList) {
			String chrom = observation.getChrom();
			if (lastChrom == null) {
				tmp = new ArrayList();
				lastChrom = chrom;
			}
			if (lastChrom != null && !lastChrom.equals(chrom)) {
				observations.add(tmp);
				tmp = new ArrayList();
				lastChrom = chrom;
			}
			tmp.add(observation);
		}
		observations.add(tmp);
		return observations;
	}

	private List<Observation> getObservationListOfChrom(List<Observation> observationList, String chrom) {
		List<List<Observation>> observationListByChrom = getObservationListByChrom(observationList);
		for (int i = 0; i < observationListByChrom.size(); i++) {
			if (observationListByChrom.get(i).get(0).getChrom().equals(chrom)) {
				return observationListByChrom.get(i);
			}
		}
		return null;
	}

	private List<Interval> getIntervalListOfChrom(List<Interval> intervalList, String chrom) {
		List<Interval> newIntervalList = new ArrayList();
		for (Interval interval : intervalList) {
			if (interval.getChrom().equals(chrom)) {
				newIntervalList.add(interval);
			}
		}
		return newIntervalList;
	}

	private List<List<Observation>> excludeRegions(List<List<Observation>> observationList,
			List<Interval> intervalList) {
		List<List<Observation>> newObservationList = new ArrayList();
		for (List<Observation> list : observationList) {
			List<Observation> currentList = new ArrayList();
			String chrom = list.get(0).getChrom();
			List<Interval> excludeList = getIntervalListOfChrom(intervalList, chrom);
			int i = 0;
			int j = 0;
			for (Interval interval : excludeList) {
				while (i < list.size() && list.get(i).getEnd() < interval.getStart()) {
					currentList.add(list.get(i));
					i++;
				}
				while (i < list.size() && list.get(i).getStart() >= interval.getStart()
						&& list.get(i).getStart() <= interval.getEnd()) {
					i++;
				}
			}
			while (i < list.size()) {
				currentList.add(list.get(i));
				i++;
			}
			newObservationList.add(currentList);
		}

		return newObservationList;
	}

	public double[] getPI() {
		int numStates = 9;
		double[] PI = new double[numStates];
		for (int i = 0; i < numStates; i++) {
			PI[i] = Math.log(1.0 / 9.0);
		}
		return PI;
	}

	private double testOverdispersion(List<Observation> observationList) throws IOException {
		int size = observationList.size();
		List<Integer> refAllele = new ArrayList();
		List<Integer> altAllele = new ArrayList();
		for (int i = 0; i < size; i++) {
			Observation observation = observationList.get(i);
			List<int[]> snvList = observation.getNormalAlleleCountList();
			for (int[] allele : snvList) {
				refAllele.add(allele[0]);
				altAllele.add(allele[1]);
			}
		}
		int[] ni = new int[refAllele.size()];
		int[] xi = new int[refAllele.size()];
		double[] pi = new double[refAllele.size()];
		for (int i = 0; i < refAllele.size(); i++) {
			ni[i] = refAllele.get(i).intValue() + altAllele.get(i).intValue();
			xi[i] = refAllele.get(i).intValue();
			pi[i] = (double) xi[i] / (double) ni[i];
		}
		Process pid = Runtime.getRuntime().exec("which Rscript");
		BufferedReader runTimeReader = new BufferedReader(new InputStreamReader(pid.getInputStream()));
		String R_HOME = runTimeReader.readLine().trim();
		String[] args = new String[] { this.afFile };
		double[] outliers = new double[] { outlier / 2, 1 - outlier / 2 };
		RCaller caller = new RCaller();
		RCode code = new RCode();
		caller.setRscriptExecutable(R_HOME);
		caller.setRCode(code);
		code.addDoubleArray("outliers", outliers);
		code.addIntArray("ni", ni);
		code.addIntArray("xi", xi);
		code.addDoubleArray("pi", pi);
		code.addRCode("library(qcc)");
		code.addRCode("data<-data.frame(ni,xi,pi)");
		code.addRCode("remove_outliers <- function(data, na.rm = TRUE, ...){");
		code.addRCode("range <- quantile(data[,3], prob=c(outliers[1], outliers[2]),na.rm =TRUE)");
		code.addRCode("new <- data");
		code.addRCode("new[which(new[,3] >= range[1] & new[,3] <= range[2]),]");
		code.addRCode("}");
		code.addRCode("range <-quantile(data[,3], prob=c(outliers[1], outliers[2]),na.rm =TRUE)");
		code.addRCode("data<-remove_outliers(data)");
		code.addRCode("if(nrow(data)>1000){");
		code.addRCode("data<-data[sample(nrow(data),1000,replace=F),]");
		code.addRCode("}");
		code.addRCode("ni<-data[,1]");
		code.addRCode("xi<-data[,2]");
		code.addRCode("test<-qcc.overdispersion.test(xi,ni)");
		code.addRCode("pvalue<-test[3]");
		code.addRCode("result<-list(pvalue=pvalue)");
		caller.runAndReturnResult("result");
		double[] pvalue = caller.getParser().getAsDoubleArray("pvalue");
		return pvalue[0];
	}

	private double[] getDispersion(List<Observation> observationList) throws IOException {
		int size = observationList.size();
		List<Integer> refAllele = new ArrayList();
		List<Integer> altAllele = new ArrayList();
		for (int i = 0; i < size; i++) {
			Observation observation = observationList.get(i);
			List<int[]> snvList = observation.getNormalAlleleCountList();
			for (int[] allele : snvList) {
				refAllele.add(allele[0]);
				altAllele.add(allele[1]);
			}
		}
		int[] ni = new int[refAllele.size()];
		int[] xi = new int[refAllele.size()];
		double[] pi = new double[refAllele.size()];

		for (int i = 0; i < refAllele.size(); i++) {
			ni[i] = refAllele.get(i).intValue() + altAllele.get(i).intValue();
			xi[i] = refAllele.get(i).intValue();
			pi[i] = (double) xi[i] / (double) ni[i];
		}

		Process pid = Runtime.getRuntime().exec("which Rscript");
		BufferedReader runTimeReader = new BufferedReader(new InputStreamReader(pid.getInputStream()));
		String R_HOME = runTimeReader.readLine().trim();
		double[] outliers = new double[] { outlier / 2, 1 - outlier / 2 };
		RCaller caller = new RCaller();
		RCode code = new RCode();
		caller.setRscriptExecutable(R_HOME);
		caller.setRCode(code);
		code.addDoubleArray("outliers", outliers);
		code.addIntArray("ni", ni);
		code.addIntArray("xi", xi);
		code.addDoubleArray("pi", pi);
		code.addRCode("library(qcc)");
		code.addRCode("library(VGAM)");
		code.addRCode("data<-data.frame(ni,xi,pi)");
		code.addRCode("remove_outliers <- function(data, na.rm = TRUE, ...){");
		code.addRCode("range <- quantile(data[,3], prob=c(outliers[1], outliers[2]),na.rm =TRUE)");
		code.addRCode("new <- data");
		code.addRCode("new[which(new[,3] >= range[1] & new[,3] <= range[2]),]");
		code.addRCode("}");
		code.addRCode("range<-quantile(data[,3], prob=c(outliers[1], outliers[2]),na.rm =TRUE)");
		code.addRCode("data<-remove_outliers(data)");
		code.addRCode("if(nrow(data)>1000){");
		code.addRCode("data<-data[sample(nrow(data),1000,replace=F),]");
		code.addRCode("}");
		code.addRCode("ni<-data[,1]");
		code.addRCode("xi<-data[,2]");
		code.addRCode("fit<-vglm(cbind(ni-xi, xi) ~ 1, betabinomial, trace = TRUE)");
		code.addRCode("parameters<-c(range,Coef(fit))");
		code.addRCode("result<-list(parameters=parameters)");
		caller.runAndReturnResult("result");
		double[] parameters = caller.getParser().getAsDoubleArray("parameters");
		return parameters;
	}

	private NBModel[] getNBModels(List<Observation> observationList) {
		int size = observationList.size();
		int[] normalDepth = new int[size];
		int[] tumorDepth = new int[size];
		int[] gc = new int[size];
		double[] mappability = new double[size];
		for (int i = 0; i < observationList.size(); i++) {
			Observation observation = observationList.get(i);
			normalDepth[i] = observation.getRD()[0];
			tumorDepth[i] = observation.getRD()[1];
			gc[i] = observation.getGC();
			mappability[i] = observation.getMappability();
		}
		try {
			NBModel[] nbModels = new NBModel[2];
			Process pid = Runtime.getRuntime().exec("which Rscript");
			BufferedReader runTimeReader = new BufferedReader(new InputStreamReader(pid.getInputStream()));
			String R_HOME = runTimeReader.readLine().trim();
			runTimeReader.close();
			double[] outliers = new double[] { outlier / 2, 1 - outlier / 2 };
			int[] binsize = new int[] { MIN_GC_BIN_SIZE };
			for (int k = 0; k < 2; k++) {
				int[] sample = new int[] { k + 1 };
				RCaller caller = new RCaller();
				RCode code = new RCode();
				caller.setRscriptExecutable(R_HOME);
				caller.setRCode(code);
				code.addIntArray("sample", sample);
				code.addIntArray("normaldepth", normalDepth);
				code.addIntArray("tumordepth", tumorDepth);
				code.addIntArray("gc", gc);
				code.addDoubleArray("mappability", mappability);
				code.addDoubleArray("outliers", outliers);
				code.addIntArray("size", binsize);
				code.addRCode("library(MASS)");
				code.addRCode("data<-data.frame(normaldepth,tumordepth,gc,mappability)");
				code.addRCode("data<-split(data[,c(sample[1],4)],data$gc)");
				code.addRCode("theta<-rep(0,101)");
				code.addRCode("coef1<-rep(0,101)");
				code.addRCode("coef2<-rep(0,101)");
				code.addRCode("remove_outliers <- function(data, na.rm = TRUE, ...) {");
				code.addRCode("range <- quantile(data[,1], prob=c(outliers[1], outliers[2]),na.rm = TRUE)");
				code.addRCode("new <- data");
				code.addRCode("new[which(new[,1] >= range[1] &  new[,1] <= range[2]),]");
				code.addRCode("}");
				code.addRCode("for(i in 1:length(data)){");
				code.addRCode("bin<-as.data.frame(unname(data[i]))");
				code.addRCode("if(nrow(bin)>size[1]){");
				code.addRCode("bin<-remove_outliers(bin)");
				code.addRCode("index<-sample(1:nrow(bin), size = min(nrow(bin),size[1]))");
				code.addRCode("count<-bin[index,1]");
				code.addRCode("mappability<-bin[index,2]");
				code.addRCode("fit<-glm.nb(count~mappability)");
				code.addRCode("theta[as.numeric(names(data[i]))+1]<-fit$theta");
				code.addRCode("coef1[as.numeric(names(data[i]))+1]<-unname(fit$coef[1])");
				code.addRCode("coef2[as.numeric(names(data[i]))+1]<-unname(fit$coef[2])");
				code.addRCode("}");
				code.addRCode("}");
				code.addRCode("result<-list(theta=theta,coef1=coef1, coef2=coef2)");
				caller.runAndReturnResult("result");
				double[] theta = caller.getParser().getAsDoubleArray("theta");
				double[] coef1 = caller.getParser().getAsDoubleArray("coef1");
				double[] coef2 = caller.getParser().getAsDoubleArray("coef2");
				nbModels[k] = new NBModel(theta, coef1, coef2);
			}
			return nbModels;
		} catch (IOException ex) {
			ex.printStackTrace();
			return null;
		}
	}

	private boolean checkInstalledPackages(String R_HOME) {
		boolean tag = true;
		RCaller caller = new RCaller();
		caller.setRscriptExecutable(R_HOME);
		caller.cleanRCode();
		RCode code = new RCode();
		String[] packages = new String[] { "Runiversal", "MASS", "VGAM", "qcc" };
		code.addStringArray("packages", packages);
		code.addRCode("label<-c(1, 1, 1, 1)");
		code.addRCode("for(i in 1:4){");
		code.addRCode("if(!require(packages[i], character.only=TRUE)){");
		code.addRCode("label[i]=0");
		code.addRCode("}");
		code.addRCode("}");
		caller.setRCode(code);
		caller.runAndReturnResult("label");
		int[] label = caller.getParser().getAsIntArray("label");
		for (int i = 0; i < packages.length; i++) {
			if (label[i] == 0) {
				logger.error(packages[i] + " is not installed in R environment!");
				tag = false;
			}
		}
		return tag;
	}

	private CNVTYPE parseCNVType(String typeString) {
		if (typeString.equals("GERMLINE_DEL")) {
			return CNVTYPE.GERMLINE_DEL;
		}
		if (typeString.equals("GERMLINE_DUP")) {
			return CNVTYPE.GERMLINE_DUP;
		}
		if (typeString.equals("SOMATIC_GAIN")) {
			return CNVTYPE.SOMATIC_GAIN;
		}
		if (typeString.equals("SOMATIC_LOSS")) {
			return CNVTYPE.SOMATIC_LOSS;
		}
		return null;
	}

}