package cn.edu.hit.tumorcnv;
/**
*
* @author Yongzhuang Liu
*/

import java.util.Comparator;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

public class ComparableSamLocusIterator extends PeekableIterator<SamLocusIterator.LocusInfo>
		implements Comparable<ComparableSamLocusIterator> {

	private String sampleID = "";

	private Comparator<SamLocusIterator.LocusInfo> comparator;

	public ComparableSamLocusIterator(final SamLocusIterator samLocusIterator, String sampleID,
			final Comparator<SamLocusIterator.LocusInfo> comparator) {
		super((java.util.Iterator<SamLocusIterator.LocusInfo>) samLocusIterator);
		this.comparator = comparator;
		this.sampleID = sampleID;
	}

	public String getSampleID() {
		return sampleID;
	}

	@Override
	public int compareTo(ComparableSamLocusIterator that) {

		final LocusInfo thisRecord = this.peek();
		final LocusInfo thatRecord = that.peek();
		return this.comparator.compare(thisRecord, thatRecord);
	}

	@Override
	public boolean equals(final Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}

		return compareTo((ComparableSamLocusIterator) o) == 0;
	}

}
