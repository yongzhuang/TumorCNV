package cn.edu.hit.tumorcnv;

/**
*
* @author Yongzhuang Liu
*/

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import java.util.Collection;
import java.util.Comparator;
import java.util.PriorityQueue;

public class MultiSamLocusIterator {

	private final PriorityQueue<ComparableSamLocusIterator> pq;

	public MultiSamLocusIterator(Collection<SamReader> readers, IntervalList intervalList) {
		this.pq = new PriorityQueue<ComparableSamLocusIterator>(readers.size());
		for (SamReader reader : readers) {
			String sampleID = "";
			for (int i = 0; i < reader.getFileHeader().getReadGroups().size(); i++) {
				SAMReadGroupRecord temp = reader.getFileHeader().getReadGroups().get(i);
				sampleID = temp.getSample();
			}
			addIfNotEmpty(new ComparableSamLocusIterator(new SamLocusIterator(reader, intervalList), sampleID,
					new Comparator<SamLocusIterator.LocusInfo>() {
						@Override
						public int compare(SamLocusIterator.LocusInfo l1, SamLocusIterator.LocusInfo l2) {
							return l1.toString().compareTo(l2.toString());
						}
					}));
		}
	}

	private void addIfNotEmpty(final ComparableSamLocusIterator iterator) {
		if (iterator.hasNext()) {
			pq.offer(iterator);
		} else {
			iterator.close();
		}
	}

	public boolean hasNext() {
		return !this.pq.isEmpty();
	}

	public LocusInfo next() {
		ComparableSamLocusIterator iterator = this.pq.poll();
		LocusInfo record = iterator.next();
		addIfNotEmpty(iterator);
		return record;
	}

	public String getCurrentIteratorSampleID() {
		ComparableSamLocusIterator iterator = this.pq.peek();
		String sampleID = iterator.getSampleID();
		return sampleID;
	}

	public void close() {
		for (ComparableSamLocusIterator iterator : pq) {
			iterator.close();
		}
	}
}
