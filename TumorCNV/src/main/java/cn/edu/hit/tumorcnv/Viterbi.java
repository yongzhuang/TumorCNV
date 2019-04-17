package cn.edu.hit.tumorcnv;

import java.util.HashMap;
import java.util.List;

/**
 * @author Yongzhuang Liu
 *
 */
public class Viterbi {

	private List<Observation> observations;
	private Emission emission;
	private Transition transition;
	private double[] PI;

	private static final int NUM_STATES = 9;
	private HashMap<Integer, int[]> STATE_MAP;;

	public Viterbi(List<Observation> observations, double[] PI, Transition transition, Emission emission) {
		this.observations = observations;
		this.PI = PI;
		this.transition = transition;
		this.emission = emission;
		init();
	}

	private void init() {
		STATE_MAP = new HashMap();
		STATE_MAP.put(0, new int[] { 0, 0 });
		STATE_MAP.put(1, new int[] { 1, 1 });
		STATE_MAP.put(2, new int[] { 2, 2 });
		STATE_MAP.put(3, new int[] { 3, 3 });
		STATE_MAP.put(4, new int[] { 4, 4 });
		STATE_MAP.put(5, new int[] { 2, 0 });
		STATE_MAP.put(6, new int[] { 2, 1 });
		STATE_MAP.put(7, new int[] { 2, 3 });
		STATE_MAP.put(8, new int[] { 2, 4 });
	}

	public int[][] viterbi() throws Exception {
		int OS = observations.size();
		int[][] path = new int[OS][NUM_STATES];
		double[][] r = new double[OS][NUM_STATES];
		for (int l = 0; l < NUM_STATES; l++) {
			int n = STATE_MAP.get(l)[0];
			int t = STATE_MAP.get(l)[1];
			double emissionProb = emission.getProbOfEmission(n, t, observations.get(0));
			r[0][l] = PI[l] + emissionProb;
			path[0][l] = l;
		}
		for (int os = 1; os < OS; os++) {
			for (int l = 0; l < NUM_STATES; l++) {
				int n = STATE_MAP.get(l)[0];
				int t = STATE_MAP.get(l)[1];
				double tmp = Double.NEGATIVE_INFINITY;
				int s = -1;
				double emissionProb = emission.getProbOfEmission(n, t, observations.get(os));
				Observation observation = observations.get(os);

				int d = observations.get(os).getStart() - observations.get(os - 1).getEnd();
				double[][] transitionProb = transition.getTransitionMatrix(d);

				List<int[]> normalAlleleList = observation.getNormalAlleleCountList();
				List<int[]> tumorAlleleList = observation.getTumorAlleleCountList();
				for (int l1 = 0; l1 < NUM_STATES; l1++) {
					double tem = r[os - 1][l1] + emissionProb + transitionProb[l1][l];
					if (tem > tmp) {
						tmp = tem;
						s = l1;
					}
				}
				r[os][l] = tmp;
				path[os][l] = s;
			}
		}
		double p = Double.NEGATIVE_INFINITY;
		int s = -1;
		for (int i = 0; i < NUM_STATES; i++) {
			if (r[r.length - 1][i] > p) {
				p = r[r.length - 1][i];
				s = i;
			}
		}
		int[] trace = new int[OS];
		trace[OS - 1] = s;
		for (int os = OS - 1; os > 0; os--) {
			trace[os - 1] = path[os][s];
			s = path[os][s];
		}
		int[][] trace2 = new int[OS][2];
		for (int i = 0; i < OS; i++) {
			int n = STATE_MAP.get(trace[i])[0];
			int t = STATE_MAP.get(trace[i])[1];
			trace2[i][0] = n;
			trace2[i][1] = t;
		}
		return trace2;
	}

	public List<Observation> getObservations() {
		return observations;
	}

	public double[] getPI() {
		return PI;
	}

	public void setPI(double[] pI) {
		PI = pI;
	}

	public Transition getTransition() {
		return transition;
	}

	public void setTransition(Transition transition) {
		this.transition = transition;
	}
}
