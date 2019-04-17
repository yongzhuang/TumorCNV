package cn.edu.hit.tumorcnv;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Yongzhuang Liu
 */
public class Transition {

	private double p;
	private int D;

	public Transition(int D, double p) {
		this.p = p;
		this.D = D;
	}

	public double[][] getTransitionMatrix(int d) {
		double deviation = 1 - Math.pow(Math.E, 0 - (double) d / (double) D);
		int numStates = 9;
		double[][] a = new double[numStates][numStates];
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				if (i == j) {
					a[i][j] = Math.log(1 - 8 * p * deviation);
				} else {
					a[i][j] = Math.log(p * deviation);
				}
			}
		}
		return a;
	}
}
