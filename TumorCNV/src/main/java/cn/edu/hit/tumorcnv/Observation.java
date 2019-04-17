package cn.edu.hit.tumorcnv;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Yongzhuang Liu
 *
 */
public class Observation {

	private String chrom;
	private int start;
	private int end;
	private int GC;
	private double mappability;


	private List<int[]> normalAlleleCountList;
	private List<int[]> tumorAlleleCountList;
	private int[] RD;

	public Observation(String chrom, int start, int end, int GC, double mappability, List<int[]> normalAlleleCountList,
			List<int[]> tumorAlleleCountList, int[] RD) {
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.GC = GC;
		this.mappability = mappability;
		this.RD = RD;
		this.normalAlleleCountList = normalAlleleCountList;
		this.tumorAlleleCountList = tumorAlleleCountList;
		if (this.normalAlleleCountList == null) {
			this.normalAlleleCountList = new ArrayList();
		}
		if (this.tumorAlleleCountList == null) {
			this.tumorAlleleCountList = new ArrayList();
		}
	}

	public void addNormalAlleleCountList(int[] alleleCount) {
		this.normalAlleleCountList.add(alleleCount);
	}

	public void addTumorAlleleCountList(int[] alleleCount) {
		this.tumorAlleleCountList.add(alleleCount);
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

	public int getGC() {
		return GC;
	}

	public double getMappability() {
		return mappability;
	}

	public List<int[]> getNormalAlleleCountList() {
		return normalAlleleCountList;
	}

	public List<int[]> getTumorAlleleCountList() {
		return tumorAlleleCountList;
	}

	public int[] getRD() {
		return RD;
	}
	
	public void setNormalAlleleCountList(List<int[]> normalAlleleCountList) {
		this.normalAlleleCountList = normalAlleleCountList;
	}

	public void setTumorAlleleCountList(List<int[]> tumorAlleleCountList) {
		this.tumorAlleleCountList = tumorAlleleCountList;
	}
}
