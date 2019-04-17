package cn.edu.hit.tumorcnv;

/**
 *
 * @author Yongzhuang Liu
 */
public class Interval {

	private String chrom;
	private int start;
	private int end;

	public Interval(String chrom, int start, int end) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
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
}
