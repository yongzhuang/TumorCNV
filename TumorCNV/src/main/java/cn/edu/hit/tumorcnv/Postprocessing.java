
package cn.edu.hit.tumorcnv;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Yongzhuang Liu
 *
 */
public class Postprocessing {

	private List<Observation> observations;
	private int[][] trace;
	private int minDistance;
	private int tau;

	public Postprocessing(List<Observation> observations, int[][] trace, int minDistance, double tau) {
		this.observations = observations;
		this.trace = trace;
		this.minDistance = minDistance;
		this.tau = (int) tau;
	}

	public List<CNVRecord> process() {
		List<CNVRecord> records = new ArrayList<CNVRecord>();
		String mergeChrom = null;
		int mergeStart = -1;
		int mergeEnd = -1;
		String mergeState = null;
		if (trace == null) {
			return records;
		}
		for (int i = 0; i < trace.length; i++) {
			if (trace[i][0] != 2 || trace[i][1] != tau) {
				String chrom = observations.get(i).getChrom();
				int start = observations.get(i).getStart();
				int end = observations.get(i).getEnd();
				String state = Integer.toString(trace[i][0]) + trace[i][1];
				if (mergeState == null) {
					mergeChrom = chrom;
					mergeStart = start;
					mergeEnd = end;
					mergeState = state;
				} else {
					if (chrom.equals(mergeChrom) && (start - mergeEnd - 1) == 0 && state.equals(mergeState)) {
						mergeEnd = end;
					} else {
						int[] states = new int[] { Integer.parseInt(mergeState.substring(0, 1)),
								Integer.parseInt(mergeState.substring(1, 2)) };
						CNVRecord record = new CNVRecord(mergeChrom, mergeStart, mergeEnd, states);
						records.add(record);
						mergeChrom = chrom;
						mergeStart = start;
						mergeEnd = end;
						mergeState = state;
					}
				}
			}
		}
		if (mergeChrom != null) {
			int[] states = new int[] { Integer.parseInt(mergeState.substring(0, 1)),
					Integer.parseInt(mergeState.substring(1, 2)) };
			records.add(new CNVRecord(mergeChrom, mergeStart, mergeEnd, states));
		}
		List<CNVRecord> cnvList = merge(records);
		List<CNVRecord> cnvList2 = updateCNVType(cnvList);
		return cnvList2;
	}

	private List<CNVRecord> merge(List<CNVRecord> records) {
		List<CNVRecord> result = new ArrayList<CNVRecord>();
		List<CNVRecord> cluster = new ArrayList<CNVRecord>();
		CNVRecord lastRecord = null;
		for (CNVRecord record : records) {
			if (lastRecord == null) {
				lastRecord = record;
				cluster.add(record);
			} else {
				int distance= record.getStart() - lastRecord.getEnd() - 1;
				if ( distance <= minDistance) {
					cluster.add(record);
					lastRecord = record;
				} else {
					result.add(mergeByCluster(cluster));
					cluster.clear();
					cluster.add(record);
					lastRecord = record;
				}
			}
		}
		if (lastRecord != null) {
			result.add(mergeByCluster(cluster));
		}
		return result;
	}

	private CNVRecord mergeByCluster(List<CNVRecord> cluster) {
		if (cluster.size() == 1) {
			return cluster.get(0);
		}
		int[][] length = new int[2][5];
		for (CNVRecord record : cluster) {
			for (int i = 0; i < 2; i++) {
				length[i][record.getStates()[i]] += record.getLength();
			}
		}
		int cnvLength = 0;
		for (int i = 0; i < cluster.size(); i++) {
			cnvLength = cnvLength + cluster.get(i).getLength();
		}
		int clusterLength = cluster.get(cluster.size() - 1).getEnd() - cluster.get(0).getStart() + 1;
		length[0][2] += clusterLength - cnvLength;
		length[1][tau] += clusterLength - cnvLength;
		int[] states = new int[2];
		for (int i = 0; i < 2; i++) {
			int index = -1;
			int max = 0;
			for (int j = 0; j < 5; j++) {
				if (length[i][j] > max) {
					max = length[i][j];
					index = j;
				}
			}
			states[i] = index;
		}
		String chrom = cluster.get(0).getChrom();
		int start = cluster.get(0).getStart();
		int end = cluster.get(cluster.size() - 1).getEnd();
		return new CNVRecord(chrom, start, end, states);
	}

	private List<CNVRecord> updateCNVType(List<CNVRecord> cnvList) {
		List<CNVRecord> newList = new ArrayList();
		for (CNVRecord cnvRecord : cnvList) {
			CNVTYPE type = null;
			int normalState = cnvRecord.getStates()[0];
			int tumorState = cnvRecord.getStates()[1];
			if (normalState == 2 && (tumorState == tau)) {
				type = CNVTYPE.REFERENCE;
			} else if (normalState == 2 && (tumorState < tau)) {
				type = CNVTYPE.SOMATIC_LOSS;
			} else if (normalState == 2 && (tumorState > tau)) {
				type = CNVTYPE.SOMATIC_GAIN;
			} else if (normalState > 2) {
				type = CNVTYPE.GERMLINE_DUP;
			} else if (normalState < 2) {
				type = CNVTYPE.GERMLINE_DEL;
			} else {
				type = CNVTYPE.NULL;
			}
			if (!type.equals(CNVTYPE.NULL) && !type.equals(CNVTYPE.REFERENCE)) {
				cnvRecord.setType(type);
				newList.add(cnvRecord);
			}
		}
		return newList;
	}
}
