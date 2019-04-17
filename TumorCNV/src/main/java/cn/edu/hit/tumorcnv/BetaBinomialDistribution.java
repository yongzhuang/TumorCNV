package cn.edu.hit.tumorcnv;

import org.apache.commons.math3.special.Gamma;

/**
 *
 * @author Yongzhuang Liu
 */
public class BetaBinomialDistribution {

	private double alpha;
	private double beta;

	public BetaBinomialDistribution(double pi, double theta) {
		this.alpha = pi / theta;
		this.beta = (1 - pi) / theta;
	}

	public double getLogDensity(int n, int k) {
		return (Gamma.logGamma(k + alpha) + Gamma.logGamma(n - k + beta) + Gamma.logGamma(alpha + beta)
				+ Gamma.logGamma(n + 1))
				- (Gamma.logGamma(alpha + beta + n) + Gamma.logGamma(alpha) + Gamma.logGamma(beta)
						+ Gamma.logGamma(k + 1) + Gamma.logGamma(n - k + 1));
	}

}
