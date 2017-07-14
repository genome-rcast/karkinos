/*
Copyright Hiroki Ueda

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package jp.ac.utokyo.karkinos.ploidy;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class AFCounter {

	public SummaryStatistics getSs() {
		return ss;
	}

	SummaryStatistics ss = new SummaryStatistics();
	float[] counter = new float[100];

	public void setValue(double d) {

		//
		ss.addValue(d);
		int idx = (int) Math.round(d * 100);
		if (idx < 0) {
			idx = 0;
		}
		if (idx > 99) {
			idx = 99;
		}
		counter[idx] = counter[idx] + 1;

	}

	public float[] getPeakDistance(double sd) {

		//
		int start = 50;
		int end = 100;

		List<Double> vals = new ArrayList<Double>();
		for (int n = start; n <= end; n++) {
			//
			float[] theoreticalDist = getThereticalDist(n, sd);
			double cp = crossProduct(theoreticalDist, counter);
			vals.add(cp);
		}

		int maxn = getMax(counter);

		double maxcp = 0;
		int maxratio = 0;
		int numofpeak = 0;

		int idx = 1;
		for (int n = start + 1; n < end - 1; n++) {
			//

			double val = vals.get(idx);
			double valb4 = vals.get(idx - 1);
			if (val == 0 || valb4 == 0) {
				idx++;
				continue;
			}
			if (val > maxcp) {
				maxcp = val;
				maxratio = n;

				double valafter = vals.get(idx + 1);
				boolean peak = (val > valb4) && (val > valafter);
				if (peak) {

					numofpeak++;

				}
			}
			idx++;
		}
		//

		if (Math.abs(maxn - 50) <= 4) {
			boolean centered = true;
			try {
				centered = counter[maxn] > (counter[maxratio] * 2);
			} catch (Exception ex) {

			}
			if (centered) {
				return new float[] { (float) ((maxn) * 0.01), 0 };
			}
		}

		//
		return new float[] { (float) ((maxratio) * 0.01), numofpeak };

	}

	private int getMax(float[] counter2) {
		float maxcnt = 0;
		int maxidx = 0;
		for (int n = 0; n < counter.length; n++) {
			if (counter[n] > maxcnt) {
				maxcnt = counter[n];
				maxidx = n;
			}
		}
		return maxidx;
	}

	private static double crossProduct(float[] theoreticalDist, float[] snpDistT) {

		float[] normalizesnpDistT = norm(snpDistT);
		double d = 0;
		for (int n = 2; n < 98; n++) {
			//
			d = d + (theoreticalDist[n] * normalizesnpDistT[n]);
		}
		return d;

	}

	private static float[] norm(float[] dist) {
		float sum = 0;
		int idx = 0;
		for (float f : dist) {
			if (idx == 0) {
				idx++;
				continue;
			}
			if (idx == 99) {
				idx++;
				continue;
			}
			if (idx >= 45 && idx <= 55) {
				idx++;
				continue;
			}
			sum = sum + f;
			idx++;
		}
		float[] ndist = new float[dist.length];

		//
		int n = 0;
		for (float f : dist) {
			if (n == 0) {
				n++;
				continue;
			}
			if (n == 99) {
				n++;
				continue;
			}
			if (n >= 45 && n <= 55) {
				n++;
				continue;
			}
			ndist[n] = (float) ((double) f / (double) sum);
			n++;
		}
		return ndist;
	}

	private static float[] getThereticalDist(int ratio, double variance) {

		float x = 0;
		double sum = 0;
		float[] dist = new float[100];
		for (int n = 0; n < 100; n++) {

			x = n + 0.5f;
			float y = dist(x, variance, ratio);
			sum = sum + y;
			dist[n] = y;
		}
		int idx = 0;
		for (double d : dist) {

			dist[idx] = (float) (dist[idx] / sum);
			idx++;
		}
		return dist;
	}

	private static float dist(float x, double sd, int ratio) {

		if (sd == 0) {
			sd = 10;// to avoid error
		}

		if (ratio == 50) {
			float mean = 50;
			NormalDistribution normald = new NormalDistributionImpl(mean, sd);
			return (float) normald.density((double) x);
		} else {

			NormalDistribution normald1 = new NormalDistributionImpl(ratio, sd);
			NormalDistribution normald2 = new NormalDistributionImpl(
					50 - (ratio - 50), sd);
			return (float) (normald1.density((double) x) + normald2
					.density((double) x));

		}

	}

}
