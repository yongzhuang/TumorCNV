package cn.edu.hit.tumorcnv;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PascalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

/**
 *
 * @author Yongzhuang Liu
 */
public class Emission {

	private NBModel normalNBModel;
	private NBModel tumorNBModel;

	private double alpha;
	private double tau;
	private double rho;

	private static final double LOG_MIN_PRO = -1000;
	private static final double LOG_MAX_PRO = -Math.pow(10, -100);
	private static final int MAX_RD = 5000;
	private static final double EPSILON = 0.1;
	private static final double DISPERSION_THRESHOLD = 500;
	private static final int MIN_READ_COUNT = 1;

	public Emission(NBModel normalNBModel, NBModel tumorNBModel, double alpha, double tau, double rho) {
		this.normalNBModel = normalNBModel;
		this.tumorNBModel = tumorNBModel;
		this.alpha = alpha;
		this.tau = tau;
		this.rho = rho;
	}

	public double getProbOfEmission(int normalState, int tumorState, Observation observation) throws Exception {
		return getProbOfRDEmission(normalState, tumorState, observation)
				+ getProbOfAFEmission(normalState, tumorState, observation);
	}

	public double getProbOfRDEmission(int normalState, int tumorState, Observation observation) throws Exception {
		double prob1 = 0;
		double prob2 = 0;
		int[] rd = observation.getRD();
		int gc = observation.getGC();
		double mappability = observation.getMappability();

		double[] theta = normalNBModel.getTheta();
		double[] coef1 = normalNBModel.getCoef1();
		double[] coef2 = normalNBModel.getCoef2();

		if (theta[gc] == 0) {
			return 0;
		}

		if (rd[0] == 0) {
			if (normalState == 0) {
				prob1 = LOG_MAX_PRO;
			} else {
				prob1 = LOG_MIN_PRO;
			}
		} else if (rd[0] > MAX_RD) {
			if (normalState == 4) {
				prob1 = LOG_MAX_PRO;
			} else {
				prob1 = LOG_MIN_PRO;
			}
		} else {
			double mu = Math.exp(coef1[gc] + coef2[gc] * mappability);
			double size = theta[gc];
			double p = 0;
			if (normalState == 0) {
				size = size * EPSILON;
				if (size > DISPERSION_THRESHOLD) {
					size = DISPERSION_THRESHOLD;
				}
				p = size / (size + mu * EPSILON);
			} else {
				size = size * (0.5 * (double) normalState);
				if (size > DISPERSION_THRESHOLD) {
					size = DISPERSION_THRESHOLD;
				}
				p = size / (size + mu * 0.5 * (double) normalState);
			}
			PascalDistribution pascal = new PascalDistribution((int) Math.ceil(size), p);
			prob1 = pascal.logProbability(rd[0]);
			if (prob1 == Double.NEGATIVE_INFINITY) {
				prob1 = LOG_MIN_PRO;
			}
		}

		theta = tumorNBModel.getTheta();
		coef1 = tumorNBModel.getCoef1();
		coef2 = tumorNBModel.getCoef2();
		if (theta[gc] == 0) {
			return 0;
		}

		if (rd[1] == 0) {
			if (tumorState == 0) {
				prob2 = LOG_MAX_PRO;
			} else {
				prob2 = LOG_MIN_PRO;
			}
		} else if (rd[1] >= MAX_RD) {
			if (tumorState == 4) {
				prob2 = LOG_MAX_PRO;
			} else {
				prob2 = LOG_MIN_PRO;
			}
		} else {
			double mu = Math.exp(coef1[gc] + coef2[gc] * mappability);
			double size = theta[gc];
			double p = 0;
			if ((normalState * (1 - alpha) + tumorState * alpha) == 0) {
				size = size * EPSILON;
				if (size > DISPERSION_THRESHOLD) {
					size = DISPERSION_THRESHOLD;
				}
				p = size / (size + mu * EPSILON);
			} else {
				size = size * (normalState * (1 - alpha) + tumorState * alpha) / (2.0 * (1 - alpha) + tau * alpha);
				if (size > DISPERSION_THRESHOLD) {
					size = DISPERSION_THRESHOLD;
				}
				p = size / (size
						+ mu * (normalState * (1 - alpha) + tumorState * alpha) / (2.0 * (1 - alpha) + tau * alpha));
			}
			PascalDistribution pascal = new PascalDistribution((int) Math.ceil(size), p);
			prob2 = pascal.logProbability(rd[1]);
			if (prob2 == Double.NEGATIVE_INFINITY) {
				prob2 = LOG_MIN_PRO;
			}
		}
		double prob = prob1 + prob2;
		if (prob == Double.POSITIVE_INFINITY || prob == Double.NEGATIVE_INFINITY || Double.isNaN(prob)) {
			throw new EmissionException();
		} else
			return prob;
	}

	public double getProbOfAFEmission(int normalState, int tumorState, Observation observation) throws Exception {
		double prob1 = 0;
		double prob2 = 0;
		List<double[]> frequencyTable = new ArrayList();
		frequencyTable.add(new double[] { 0.01 });
		frequencyTable.add(new double[] { 0.5 });
		frequencyTable.add(new double[] { 0.33, 0.67 });
		frequencyTable.add(new double[] { 0.25, 0.5, 0.75 });
		List<double[]> frequencyTable2 = new ArrayList();
		frequencyTable2.add(new double[] { 0 });
		frequencyTable2.add(new double[] { 0.01, 0.99 });
		frequencyTable2.add(new double[] { 0.01, 0.5, 0.99 });
		frequencyTable2.add(new double[] { 0.01, 0.33, 0.67, 0.99 });
		frequencyTable2.add(new double[] { 0.01, 0.25, 0.5, 0.75, 0.99 });
		List<int[]> normalAlleleCountList = observation.getNormalAlleleCountList();
		List<int[]> tumorAlleleCountList = observation.getTumorAlleleCountList();
		if (normalAlleleCountList.size() == 0 || normalAlleleCountList.size() > 3) {
			return 0;
		} else if (normalState == 0) {
			return LOG_MIN_PRO;
		}
		if (normalState == 1) {
			for (int[] alleleCount : normalAlleleCountList) {
				if (rho == 0) {
					BinomialDistribution b = new BinomialDistribution(alleleCount[0] + alleleCount[1], 0.01);
					prob1 += b.logProbability(alleleCount[1]);
				} else {
					BetaBinomialDistribution bb = new BetaBinomialDistribution(0.01, rho);
					prob1 += bb.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[1]);
				}
			}
		} else if (normalState == 2) {
			for (int[] alleleCount : normalAlleleCountList) {
				if (rho == 0) {
					BinomialDistribution b = new BinomialDistribution(alleleCount[0] + alleleCount[1], 0.5);
					prob1 += b.logProbability(alleleCount[0]);
				} else {
					BetaBinomialDistribution bb = new BetaBinomialDistribution(0.5, rho);
					prob1 += bb.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[0]);
				}
			}
		} else if (normalState == 3) {
			for (int[] alleleCount : normalAlleleCountList) {
				double[] tmpProb = new double[2];
				if (rho == 0) {
					BinomialDistribution b1 = new BinomialDistribution(alleleCount[0] + alleleCount[1], 0.33);
					tmpProb[0] = b1.logProbability(alleleCount[0]);
					BinomialDistribution b2 = new BinomialDistribution(alleleCount[0] + alleleCount[1], 0.67);
					tmpProb[1] = b2.logProbability(alleleCount[0]);
				} else {
					BetaBinomialDistribution bb1 = new BetaBinomialDistribution(0.33, rho);
					tmpProb[0] = bb1.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[0]);
					BetaBinomialDistribution bb2 = new BetaBinomialDistribution(0.67, rho);
					tmpProb[1] = bb2.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[0]);
				}
				prob1 += getAverage(tmpProb);
			}
		} else if (normalState == 4) {
			for (int[] alleleCount : normalAlleleCountList) {
				double[] tmpProb = new double[3];
				if (rho == 0) {
					BinomialDistribution b1 = new BinomialDistribution(alleleCount[0] + alleleCount[1], 0.25);
					tmpProb[0] = b1.logProbability(alleleCount[0]);
					BinomialDistribution b2 = new BinomialDistribution(alleleCount[0] + alleleCount[1], 0.5);
					tmpProb[1] = b2.logProbability(alleleCount[0]);
					BinomialDistribution b3 = new BinomialDistribution(alleleCount[0] + alleleCount[1], 0.75);
					tmpProb[2] = b3.logProbability(alleleCount[0]);
				} else {
					BetaBinomialDistribution bb1 = new BetaBinomialDistribution(0.25, rho);
					tmpProb[0] = bb1.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[0]);
					BetaBinomialDistribution bb2 = new BetaBinomialDistribution(0.5, rho);
					tmpProb[1] = bb2.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[0]);
					BetaBinomialDistribution bb3 = new BetaBinomialDistribution(0.75, rho);
					tmpProb[2] = bb3.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[0]);
				}
				prob1 += getAverage(tmpProb);
			}
		}

		for (int[] alleleCount : tumorAlleleCountList) {
			double total_copy_number_state = (double) normalState * (1 - alpha) + (double) tumorState * alpha;
			if ((alleleCount[0] + alleleCount[1]) < MIN_READ_COUNT) {
				if (total_copy_number_state == 0) {
					prob2 += LOG_MAX_PRO;
				} else {
					prob2 += LOG_MIN_PRO;
				}
			} else {
				if (total_copy_number_state == 0) {
					prob2 += LOG_MIN_PRO;
				} else {
					double tmp = 0;
					double[] tmpProb = new double[frequencyTable.get(normalState - 1).length
							* frequencyTable2.get(tumorState).length];
					int index = 0;
					for (int m = 0; m < frequencyTable.get(normalState - 1).length; m++) {
						for (int n = 0; n < frequencyTable2.get(tumorState).length; n++) {
							double p1 = frequencyTable.get(normalState - 1)[m];
							double p2 = frequencyTable2.get(tumorState)[n];
							double p = ((double) normalState * p1 * (1 - alpha) + (double) tumorState * p2 * alpha)
									/ total_copy_number_state;
							if (rho == 0) {
								BinomialDistribution b = new BinomialDistribution(alleleCount[0] + alleleCount[1], p);
								tmpProb[index] = b.logProbability(alleleCount[0]);
							} else {
								BetaBinomialDistribution bb = new BetaBinomialDistribution(p, rho);
								tmpProb[index] = bb.getLogDensity(alleleCount[0] + alleleCount[1], alleleCount[0]);
							}
							index++;
						}
					}
					prob2 += getAverage(tmpProb);
				}
			}
		}
		double prob = prob1 + prob2;
		if (prob == Double.POSITIVE_INFINITY || prob == Double.NEGATIVE_INFINITY || Double.isNaN(prob)) {
			throw new EmissionException();

		} else {
			return prob;
		}
	}

	private double getAverage(double[] array) throws EmissionException {
		double prob = getLogSum(array) - Math.log(array.length);
		if (prob == Double.POSITIVE_INFINITY || prob == Double.NEGATIVE_INFINITY || Double.isNaN(prob)) {
			System.out.println("prob=" + prob);
			throw new EmissionException();

		} else
			return prob;
	}

	public double getLogSum(double[] xs) {
		if (xs.length == 1)
			return xs[0];
		double max = xs[0];
		for (int i = 0; i < xs.length; i++) {
			if (max < xs[i]) {
				max = xs[i];
			}
		}
		double sum = 0.0;
		for (int i = 0; i < xs.length; ++i)
			if (xs[i] != Double.NEGATIVE_INFINITY)
				sum += java.lang.Math.exp(xs[i] - max);
		return max + java.lang.Math.log(sum);
	}
}
