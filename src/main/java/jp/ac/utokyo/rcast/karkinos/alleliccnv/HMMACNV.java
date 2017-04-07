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
package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.DistMedian;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;

public class HMMACNV {

	public static void calc(List<List<SNVHolderPlusACnv>> plist)
			throws IOException {

		
		DistMedian low = new DistMedian(0.5,1.5,true);
		for(List<SNVHolderPlusACnv> list : plist){
			for (SNVHolderPlusACnv sc : list) {

				low.addValue(sc.getLowera().getWtval());

			}
		}
		double d = low.getDistributionMedian();
		double adjust = 0;
		if(d>0.5&&d<1){
			adjust = 1-d;
		}
		
		Hmm<ObservationReal> hmmhigh = getHMM(plist, true,0);
		Hmm<ObservationReal> hmmlow = getHMM(plist, false,adjust);

		for (List<SNVHolderPlusACnv> list : plist) {

			//
			List<ObservationReal> olisthigh = getList(list, true);
			List<ObservationReal> olistlow = getList(list, false);
			int[] hmmaryhigh = initAry(list.size());
			int[] hmmarylow = initAry(list.size());
			try {
				hmmaryhigh = hmmhigh.mostLikelyStateSequence(olisthigh);
			} catch (IllegalArgumentException ea) {

			}
			try {
				hmmarylow = hmmlow.mostLikelyStateSequence(olistlow);
			} catch (IllegalArgumentException ea) {

			}
			int cnt = 0;
			for (SNVHolderPlusACnv sc : list) {

				//
				sc.getHighera().setHmmval(hmmaryhigh[cnt] + 1);
				sc.getLowera().setHmmval(hmmarylow[cnt] + 1);
				cnt++;

			}

		}

	}

	private static int[] initAry(int i) {
		
		int[] ary = new int[i];
		for(int n=0;n<i;n++){
			ary[n]=1;
		}
		return ary;
		
	}

	private static List<ObservationReal> getList(List<SNVHolderPlusACnv> list,
			boolean high) {
		List<ObservationReal> olist = new ArrayList<ObservationReal>();
		for (SNVHolderPlusACnv wi : list) {

			double d = 0;
			if (high) {
				d = wi.getHighera().getWtval();
			} else {
				d = wi.getLowera().getWtval();
			}
			ObservationReal oi = new ObservationReal(d);
			olist.add(oi);
		}
		return olist;
	}

	private static Hmm<ObservationReal> getHMM(
			List<List<SNVHolderPlusACnv>> plist, boolean high, double adjust) {

		int[] countn = new int[10];
		Set<Double> checkReg = new HashSet<Double>();

		double baseline = getBaseLine(plist, high,adjust);
		double stepSize = getStepSize(plist, high);

		int total = 0;
		int nodesize = 0;
		for (int m : countn) {
			if (m > 0) {
				nodesize++;
			}
			total = total + m;
		}
		// asuume normal distribution for each n=1,2,3,4
		OpdfGaussianFactory factory = new OpdfGaussianFactory();
		Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(4, factory);
		int idx = 0;
		for (int m : countn) {

			double p = (double) m / (double) total;
			System.out.println(idx + "\t" + p);
			if (idx >= nodesize)
				break;
			if ((p == 0) || Double.isNaN(p)) {
				if (idx == 1) {
					p = 0.9;
				} else {
					p = 0.1;
				}
			}
			hmm.setPi(idx, p);
			idx++;
		}
		int multifuctor = 20;

		double variance = Math.max((1 - baseline), stepSize) / 3;
		// 1-0.2,
		// 0.5-0.1
		System.out.println("high=" + high + " baseline=" + baseline
				+ " stepsize=" + stepSize + " variance=" + variance);

		hmm.setOpdf(0, new OpdfGaussian((baseline-adjust), variance));
		hmm.setOpdf(1, new OpdfGaussian((1-adjust), variance));
		hmm.setOpdf(2, new OpdfGaussian((1-adjust) + stepSize, variance));
		hmm.setOpdf(3, new OpdfGaussian((1-adjust) + (2 * stepSize), variance));

		hmm.setAij(0, 1, 0.1);
		hmm.setAij(0, 0, 0.9);
		hmm.setAij(0, 2, 0);
		hmm.setAij(0, 3, 0);

		hmm.setAij(1, 0, 0.05);
		hmm.setAij(1, 1, 0.9);
		hmm.setAij(1, 2, 0.05);
		hmm.setAij(1, 3, 0);

		hmm.setAij(2, 0, 0);
		hmm.setAij(2, 1, 0.08);
		hmm.setAij(2, 2, 0.9);
		hmm.setAij(2, 3, 0.2);

		hmm.setAij(3, 0, 0);
		hmm.setAij(3, 1, 0);
		hmm.setAij(3, 2, 0.1);
		hmm.setAij(3, 3, 0.9);

		return hmm;
	}

	private static SummaryStatistics getVariance(
			List<List<SNVHolderPlusACnv>> plist, boolean high) {
		SummaryStatistics ss = new SummaryStatistics();
		for (List<SNVHolderPlusACnv> list : plist) {
			for (SNVHolderPlusACnv sc : list) {
				double d = 0;
				if (high) {
					d = sc.getHighera().getWtval();
				} else {
					d = sc.getLowera().getWtval();
				}
				if (between(d, 0.8, 1.2)) {
					ss.addValue(d);
				}
			}
		}
		return ss;
	}

	private static boolean between(double d, double s, double e) {

		//
		return (d >= s) && (d <= e);

	}

	private static double getStepSize(List<List<SNVHolderPlusACnv>> plist,
			boolean high) {

		double s = 1.3;
		double e = 2;
		double sl = getMinSDLine(plist, true, s, e);
		if ((s == sl) || (e == sl)) {
			sl = 2;
		}
		return Math.abs(1 - sl);
	}

	private static double getBaseLine(List<List<SNVHolderPlusACnv>> plist,
			boolean high, double adjust) {

		double s = 0.1;
		double e = 0.7-adjust;
		double bl = getMinSDLine(plist, false, s, e);
		if ((s == bl) || (e == bl)) {
			bl = 0.1;
		}
		return bl;
	}

	private static double getMinSDLine(List<List<SNVHolderPlusACnv>> plist,
			boolean high, double s, double e) {

		List<Double> vals = new ArrayList<Double>();
		for (List<SNVHolderPlusACnv> list : plist) {
			for (SNVHolderPlusACnv sc : list) {
				double d = 0;
				if (high) {
					d = sc.getHighera().getWtval();
				} else {
					d = sc.getLowera().getWtval();
				}
				if (between(d, s, e)) {
					vals.add(d);
				}
			}
		}
		double step = 0.01;
		double minsd = 10000;
		double minval = 0;
		for (double d = s; d <= e; d = d + step) {

			double sd = dist(vals, d);
			if (sd < minsd) {
				minsd = sd;
				minval = d;
			}
		}
		return minval;
	}

	private static double dist(List<Double> vals, double d) {
		double sum = 0;
		for (double dd : vals) {
			double dif = d - dd;
			double d2 = Math.pow(dif, 2);
			sum = sum + d2;
		}
		return sum;
	}
}
