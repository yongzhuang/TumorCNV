package cn.edu.hit.tumorcnv;

/**
 *
 * @author Yongzhuang Liu
 */

public class EmissionException extends Exception {
	public EmissionException() {
		super("Emission probability is out of bounds!");
	}
}
