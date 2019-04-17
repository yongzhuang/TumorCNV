package cn.edu.hit.tumorcnv;

/**
 *
 * @author Yongzhuang Liu
 */
public class CNVRecord implements Comparable<CNVRecord> {

	private String chrom;
	private int start;
	private int end;
	private int[] states;
	private int length;
	private double mappability = -1;
	private int[] readDepth = null;
	CNVTYPE type = null;

	public CNVRecord(String chrom, int start, int end, CNVTYPE type) {
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.type = type;
		this.length = end - start + 1;
	}

	public CNVRecord(String chrom, int start, int end, CNVTYPE type, double mappability) {
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.type = type;
		this.length = end - start + 1;
		this.mappability = mappability;
	}

	public CNVRecord(String chrom, int start, int end, int[] states, double mappability) {
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.states = states;
		this.length = end - start + 1;
		this.mappability = mappability;
	}

	public CNVRecord(String chrom, int start, int end, int[] states, CNVTYPE type, double mappability) {
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.states = states;
		this.type = type;
		this.length = end - start + 1;
		this.mappability = mappability;
	}

	public CNVRecord(String chrom, int start, int end, int[] states) {
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.states = states;
		this.length = end - start + 1;
	}

	public String getChrom() {
		return chrom;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public int[] getStates() {
		return states;
	}

	public int getLength() {
		return length;
	}

	public double getMappability() {
		return mappability;
	}

	public int compareTo(CNVRecord other) {
		if (chrom.equals(other.getChrom())) {
			return start - other.start;
		} else {
			return chrom.compareTo(other.getChrom());
		}
	}

	public int[] getReadDepth() {
		return readDepth;
	}

	public void setReadDepth(int[] readDepth) {
		this.readDepth = readDepth;
	}

	public String toString() {
		return chrom + "\t" + start + "\t" + end + "\t" + this.states[0] + "\t" + this.states[1] + "\t" + toTypeString()
				+ "\t" + mappability + "\t" + length;
	}

	public CNVTYPE getType() {
		return type;
	}

	public void setType(CNVTYPE type) {
		this.type = type;
	}

	public String toTypeString() {
		if (type == null) {
			return null;
		} else {
			if (type.equals(CNVTYPE.GERMLINE_DEL)) {
				return "GERMLINE_DEL";
			} else if (type.equals(CNVTYPE.GERMLINE_DUP)) {
				return "GERMLINE_DUP";
			} else if (type.equals(CNVTYPE.SOMATIC_GAIN)) {
				return "SOMATIC_GAIN";
			} else if (type.equals(CNVTYPE.SOMATIC_LOSS)) {
				return "SOMATIC_LOSS";
			}
		}
		return null;
	}
}
