package cn.edu.hit.tumorcnv;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.util.List;

/**
 * @author Yongzhuang Liu
 *
 */
public class SingleThreadCall implements Runnable {
	private Viterbi viterbi;
	private String outputFile;
	private int minDistance;
	private double tau;

	public SingleThreadCall(Viterbi viterbi, String outputFile, int minDistance,
			double tau) {
		this.viterbi = viterbi;
		this.outputFile = outputFile;
		this.minDistance = minDistance;
		this.tau = tau;
	}

	public void run() {
		try {
			int[][] trace = this.viterbi.viterbi();
			Postprocessing post = new Postprocessing(this.viterbi.getObservations(), trace, minDistance, tau);
			List<CNVRecord> list = post.process();
			FileOutputStream out = new FileOutputStream(new File(outputFile), true);
			FileChannel fcout = out.getChannel();
			FileLock flout = null;
			while (true) {
				try {
					flout = fcout.tryLock();
					break;
				} catch (OverlappingFileLockException e) {
					try {
						Thread.sleep(100);
					} catch (InterruptedException ex) {
						ex.printStackTrace();
					}
				}
			}
			for (CNVRecord line : list) {
				String tmp = line.toString() + "\n";
				out.write(tmp.getBytes());
			}
			flout.release();
			fcout.close();
			out.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
	}
}
