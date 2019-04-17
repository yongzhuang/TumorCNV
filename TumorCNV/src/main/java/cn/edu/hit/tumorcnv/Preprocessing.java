package cn.edu.hit.tumorcnv;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.apache.log4j.Logger;

/**
 *
 * @author Yongzhuang Liu
 */
public class Preprocessing {

	static Logger logger = Logger.getLogger(Preprocessing.class);

	private String normalVCFFile = null;
	private String normalBAMFile = null;
	private String tumorBAMFile = null;
	private String referenceSequenceFile = null;
	private String mappabilityFile = null;
	private String outputPrefix = null;

	private HashMap<String, String> sampleID = new HashMap<String, String>();

	private int minMappingQuality = 0;
	private int minBaseQuality = 0;
	private int windowSize;
	private static final int MIN_DEPTH = 10;

	public Preprocessing(String normalVCFFile, String normalBAMFile, String tumorBAMFile, String referenceSequenceFile,
			String mappabilityFile, String outputPrefix, int minMappingQuality, int minBaseQuality, int windowSize) {
		this.normalVCFFile = normalVCFFile;
		this.normalBAMFile = normalBAMFile;
		this.tumorBAMFile = tumorBAMFile;
		this.referenceSequenceFile = referenceSequenceFile;
		this.mappabilityFile = mappabilityFile;
		this.outputPrefix = outputPrefix;
		this.minMappingQuality = minMappingQuality;
		this.minBaseQuality = minBaseQuality;
		this.windowSize = windowSize;
		this.init();
	}

	public void run() {
		String outputFile = "";
		if (normalVCFFile != null) {
			logger.info("Getting Allele Frequency Information ...");
			outputFile = this.outputPrefix + ".AF";
			try {
				this.getAF(outputFile);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		logger.info("Getting Read Depth ...");
		outputFile = this.outputPrefix + ".RD";
		try {
			this.getRD(outputFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void init() {
		SamReader normalSAMFileReader = SamReaderFactory.makeDefault().open(new File(this.normalBAMFile));
		for (int i = 0; i < normalSAMFileReader.getFileHeader().getReadGroups().size(); i++) {
			SAMReadGroupRecord temp = normalSAMFileReader.getFileHeader().getReadGroups().get(i);
			if (this.sampleID.values().contains(temp.getSample())) {
				continue;
			} else {
				this.sampleID.put("normal", temp.getSample());
			}
		}
		SamReader tumorSAMFileReader = SamReaderFactory.makeDefault().open(new File(this.tumorBAMFile));
		for (int i = 0; i < tumorSAMFileReader.getFileHeader().getReadGroups().size(); i++) {
			SAMReadGroupRecord temp = tumorSAMFileReader.getFileHeader().getReadGroups().get(i);
			if (this.sampleID.values().contains(temp.getSample())) {
				continue;
			} else {
				this.sampleID.put("tumor", temp.getSample());
			}
		}
	}

	public void getRD(String outputFile) {
		try {
			List<SamReader> readers = new ArrayList<SamReader>();
			SamReader normalSAMFileReader = SamReaderFactory.makeDefault().open(new File(this.normalBAMFile));
			SamReader tumorSAMFileReader = SamReaderFactory.makeDefault().open(new File(this.tumorBAMFile));
			readers.add(normalSAMFileReader);
			readers.add(tumorSAMFileReader);
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
			bw.write("CHROM\tSTART\tEND\tRD_NORMAL\tRD_TUMOR\tGC\tMappability\n");
			ReferenceSequenceFileWalker walker = new ReferenceSequenceFileWalker(new File(this.referenceSequenceFile));
			int[] RD = new int[readers.size()];
			HashMap<String, Integer> indexMap = new HashMap<String, Integer>();
			int index = 0;
			for (SamReader tempReader : readers) {
				for (int i = 0; i < tempReader.getFileHeader().getReadGroups().size(); i++) {
					SAMReadGroupRecord temp = tempReader.getFileHeader().getReadGroups().get(i);
					if (indexMap.keySet().contains(temp.getSample())) {
						continue;
					} else {
						indexMap.put(temp.getSample(), index);
						index++;
					}
				}
			}
			MultiSamRecordIterator iterators = new MultiSamRecordIterator((Collection<SamReader>) readers);
			int bin = -1;
			String chrom = null;
			byte[] bases = null;
			double[] mappabilities = null;
			while (iterators.hasNext()) {
				SAMRecord nextRead = iterators.next();
				if (nextRead.getNotPrimaryAlignmentFlag() || nextRead.getReadUnmappedFlag()
						|| nextRead.getDuplicateReadFlag() || nextRead.getInferredInsertSize() == 0
						|| nextRead.getMappingQuality() < this.minMappingQuality) {
					continue;
				}
				String id = nextRead.getReadGroup().getSample();
				String referenceName = nextRead.getReferenceName();
				int referenceIndex = nextRead.getReferenceIndex();
				Pattern pattern = Pattern.compile("^(chr)?([1-9]|1[0-9]|2[0-2]|[X|Y])$");
				Matcher matcher = pattern.matcher(referenceName);
				if (!matcher.matches()) {
					break;
				}
				int start = nextRead.getAlignmentStart();
				int currentBIN = (start - 1) / windowSize;

				if (chrom == null || !referenceName.equals(chrom)) {
					bin = currentBIN;
					RD = new int[readers.size()];
					chrom = referenceName;
					bases = walker.get(referenceIndex).getBases();
					int length = bases.length;
					mappabilities = this.getMappability(this.mappabilityFile, chrom, length, this.windowSize);
				}

				int GC;
				while (currentBIN != bin) {
					GC = (new GCCalculation(bases, bin * windowSize + 1, (bin + 1) * windowSize + 1)).getGCContent();
					double mappability = mappabilities[bin];
					bw.write(buildRecord(chrom, bin, windowSize, RD, GC, mappability) + "\n");
					bin = bin + 1;
					RD = new int[readers.size()];
				}
				RD[indexMap.get(id)]++;
			}
			bw.close();
			walker.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
	}

	private String buildRecord(String chrom, int bin, int window, int[] RD, int GC, double mappability) {
		String result = chrom + "\t" + (bin * window + 1) + "\t" + ((bin + 1) * window);
		for (int i = 0; i < RD.length; i++) {
			result += "\t" + RD[i];
		}
		result += "\t" + GC + "\t" + mappability;
		return result;
	}

	private double[] getMappability(String mappabilityFilePath, String chrom, int length, int windowSize) {
		if (!chrom.startsWith("chr")) {
			chrom = "chr" + chrom;
		}
		try {
			BBFileReader reader = new BBFileReader(mappabilityFilePath);
			BigWigIterator iterator = reader.getBigWigIterator(chrom, 1, chrom, length, true);
			double[] result = new double[length / windowSize];
			int bin = 0;
			double tmp = 0;
			while (iterator.hasNext()) {
				WigItem item = iterator.next();
				int start = item.getStartBase();
				int end = item.getEndBase();
				double value = item.getWigValue();
				for (int i = start + 1; i <= end; i++) {
					int currentBin = (i - 1) / windowSize;
					if (currentBin != bin) {
						if (tmp > 0) {
							result[bin] = tmp / windowSize;
							tmp = 0;
						}
						bin = currentBin;
					}
					tmp += value;
				}
			}
			reader.close();
			return result;
		} catch (IOException ex) {
			ex.printStackTrace();
			return null;
		}
	}

	private void getAF(String outputFile) throws IOException {
		Pattern pattern = Pattern.compile("^(chr)?([1-9]|1[0-9]|2[0-2]|[X|Y])$");
		VCFFileReader vcf = new VCFFileReader(new File(this.normalVCFFile));
		SamReader normalSAMFileReader = SamReaderFactory.makeDefault().open(new File(this.normalBAMFile));
		SamReader tumorSAMFileReader = SamReaderFactory.makeDefault().open(new File(this.tumorBAMFile));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		bw.write("CHROM\tPOSITION\tREF_COUNT_NORMAL\tALT_COUNT_NORMAL\tREF_COUNT_TUMOR\tALT_COUNT_TUMOR\n");
		IntervalList normalIntervalList = new IntervalList(normalSAMFileReader.getFileHeader());
		IntervalList tumorIntervalList = new IntervalList(tumorSAMFileReader.getFileHeader());
		CloseableIterator<VariantContext> iterator = vcf.iterator();
		Map<String, String> intervalRefAltMap = new HashMap<String, String>();
		while (iterator.hasNext()) {
			VariantContext variant = iterator.next();
			Genotype genotype = variant.getGenotype(0);
			if (variant.isSNP() && genotype.isHet() && variant.getAlternateAlleles().size() == 1
					&& variant.isBiallelic() && genotype.getDP() >= MIN_DEPTH) {
				String chrom = variant.getChr();
				int start = variant.getStart();
				int end = variant.getEnd();
				String refBase = variant.getReference().getBaseString();
				String altBase = variant.getAlternateAlleles().get(0).getBaseString();
				Matcher matcher = pattern.matcher(chrom);
				if (!matcher.matches()) {
					break;
				}
				intervalRefAltMap.put(chrom + ":" + start, refBase + "|" + altBase);
				normalIntervalList.add(new Interval(chrom, start, end));
				tumorIntervalList.add(new Interval(chrom, start, end));
			}
		}
		vcf.close();
		SamLocusIterator normalSamLocusIterator = new SamLocusIterator(normalSAMFileReader, normalIntervalList);
		SamLocusIterator tumorSamLocusIterator = new SamLocusIterator(tumorSAMFileReader, tumorIntervalList);
		for (int i = 0; i < normalIntervalList.size(); i++) {
			int[] count1 = new int[2];
			if (!normalSamLocusIterator.hasNext()) {
				break;
			}
			SamLocusIterator.LocusInfo normalLocusInfo = normalSamLocusIterator.next();
			List<SamLocusIterator.RecordAndOffset> recordAndOffsetList = normalLocusInfo.getRecordAndPositions();
			for (SamLocusIterator.RecordAndOffset recordAndOffset : recordAndOffsetList) {
				byte[] refAltBaseByte = intervalRefAltMap.get(normalLocusInfo.toString()).getBytes();
				byte base = recordAndOffset.getReadBase();
				if (base == refAltBaseByte[0] && recordAndOffset.getBaseQuality() > this.minBaseQuality
						&& recordAndOffset.getRecord().getMappingQuality() > this.minMappingQuality) {
					count1[0]++;
				}
				if (base == refAltBaseByte[2] && recordAndOffset.getBaseQuality() > this.minBaseQuality
						&& recordAndOffset.getRecord().getMappingQuality() > this.minMappingQuality) {
					count1[1]++;
				}
			}
			int[] count2 = new int[2];
			if (!tumorSamLocusIterator.hasNext()) {
				break;
			}
			SamLocusIterator.LocusInfo tumorLocusInfo = tumorSamLocusIterator.next();
			List<SamLocusIterator.RecordAndOffset> tumorRecordAndOffsetList = tumorLocusInfo.getRecordAndPositions();
			for (SamLocusIterator.RecordAndOffset recordAndOffset : tumorRecordAndOffsetList) {
				byte[] refAltBaseByte = intervalRefAltMap.get(tumorLocusInfo.toString()).getBytes();
				byte base = recordAndOffset.getReadBase();
				if (base == refAltBaseByte[0] && recordAndOffset.getBaseQuality() > this.minBaseQuality
						&& recordAndOffset.getRecord().getMappingQuality() > this.minMappingQuality) {
					count2[0]++;
				}
				if (base == refAltBaseByte[2] && recordAndOffset.getBaseQuality() > this.minBaseQuality
						&& recordAndOffset.getRecord().getMappingQuality() > this.minMappingQuality) {
					count2[1]++;
				}
			}
			bw.write(normalLocusInfo.toString().replaceFirst(":", "\t") + "\t" + count1[0] + "\t" + count1[1] + "\t"
					+ count2[0] + "\t" + count2[1] + "\n");
		}
		bw.close();
	}
}
