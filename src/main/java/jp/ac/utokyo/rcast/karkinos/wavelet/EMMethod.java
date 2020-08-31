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
package jp.ac.utokyo.rcast.karkinos.wavelet;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.cntwavelet.CTWaveletBean;
import jp.ac.utokyo.rcast.karkinos.cntwavelet.GaussianWavelet;
import jp.ac.utokyo.rcast.karkinos.cntwavelet.MexicanHatWavelet;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class EMMethod {

	public static PeaksInfo calc(DataSet dataset, double lohestimate)
			throws IOException {

		PeaksInfo pi = new PeaksInfo();
		// take over all bin distribution
		int[] bininit = new int[4000];
		int[] bin = new int[4000];

		List<Peak> peaklist = new ArrayList<Peak>();

		SummaryStatistics ss = new SummaryStatistics();

		List<List<WaveletIF>> plist = dataset.getCapInterval();
		long total = 0;
		for (List<WaveletIF> dlist : plist) {

			for (int n = 0; n < dlist.size(); n++) {

				total++;
				WaveletIF wif = dlist.get(n);
				double val = wif.getDenioseValue();
				//System.out.println("denoise " + val);
				if(Double.isNaN(val))val=1.0d;
				int idx = getIdx(val);
				bininit[idx] = bininit[idx] + 1;
				if (Math.abs(1 - val) < 0.25) {
					ss.addValue(val);
				}

			}

		}
		// remove unexcepted leap in data
		for (int n = 0; n < 3998; n++) {

			int v0 = bininit[n];
			int v1 = bininit[n + 1];
			int v2 = bininit[n + 2];
			int ave = (v0 + v2) / 2;
			boolean v1isnise = v1 > (ave * 3);
			if (v1isnise) {
				bin[n + 1] = ave;
			} else {
				bin[n + 1] = v1;
			}

		}

		// take moving average
		double[] binma = new double[4000];
		for (int n = 0; n < 4000; n++) {

			double ma = 0;
			ma = ave(bin, n - 5, n + 5);
			binma[n] = ma;

		}
		// take continuous wavelet transform
		// using mexan hat function
		// from 0.2 to 3.8
		// what is
		// double[] binma2 = new double[4000];
		// for (int n = 0; n < 4000; n++) {
		// binma2[n] = 0;
		// }
		// CTWaveletBean bean = MexicanHatWavelet.getWaveletTransform(binma);
		CTWaveletBean bean = GaussianWavelet.getWaveletTransform(binma);
		double[] binma2 = bean.getData();
		// for (int n = 200; n < 3800; n++) {
		//
		// binma2[n] = getConvolve(binma,n);
		//
		// }

		// find peaks
		double peaksum = 0;
		double prevp = 0;
		for (int n = 10; n < 3989; n++) {

			// peak

			//
			double b4b = binma2[n - 2];
			double b4 = binma2[n - 1];
			double v = binma2[n];
			double next = binma2[n + 1];
			double nextnext = binma2[n + 2];

			double x = n * 0.001;
			double y = v;

			// System.out.println(x+"\t"+y);
			boolean peak = isPeak(x, b4b, b4, v, next, nextnext);

			if (v > 0.005) {

				// if ((v >= b4) && (v >= next)) {
				//
				// if ((b4b <= b4) && (next >= nextnext)) {
				if (peak) {

					boolean toonear = Math.abs(x - prevp) < (bean.getSd() * 3);
					if (toonear) {

						if (peaklist.size() > 0) {
							Peak lastp = peaklist.get(peaklist.size() - 1);
							lastp.setXYifYBigger(x, y);
							prevp = x;

						}
						continue;
					}

					Peak p = new Peak(x, y);
					p.setU(x);
					p.setV(bean.getMostfittingvariance());
					peaklist.add(p);
					peaksum = peaksum + y;
					prevp = x;
				}
				// }

			}

		}

		for (Peak peak : peaklist) {
			double r = peak.getY() / peaksum;
			peak.setR(r);
		}

		// /
		List<EMval> list = new ArrayList<EMval>();
		for (List<WaveletIF> dlist : plist) {

			for (int n = 0; n < dlist.size(); n++) {

				WaveletIF wif = dlist.get(n);
				EMval emval = new EMval();
				emval.setValue(wif.getValue());
				list.add(emval);
			}

		}

		pi.setSignalcount(bin);
		pi.setPeaksignals(binma2);
		pi.setPeaklist(peaklist);
		pi.setMa(binma);

		pi.setPeaklist(peaklist);

		double psum = 0;
		for (Peak p : peaklist) {
			psum = psum + p.getR();

		}
		int pidx = 0;
		double maxr = 0;
		double minu = 0;
		for (Peak p : peaklist) {
			double r = p.getR() / psum;
			if (r > maxr) {
				maxr = r;
			}
			
			if((minu==0)||(minu > p.getU())){
				minu = p.getU();
			}
			
			p.setR(r);
			p.setPeakidx(pidx);
			pidx++;
			System.out.println(p.getR() + "\t" + p.getU() + "\t" + p.getV());
		}

		return pi;

	}

	private static boolean isPeak(double x, double b4b, double b4, double v,
			double next, double nextnext) {

		double diff1 = b4 - b4b;
		double diff2 = v - b4;
		double diff3 = next - v;
		double diff4 = nextnext - next;
		// first derivative
		if ((diff1 > 0) && (diff2 > 0)) {

			//
			if ((diff3 < 0) && (diff4 < 0)) {
				return true;
			}

		}
		// second derivative
		double diffd2 = diff2 - diff1;
		double diffd3 = diff3 - diff2;
		double diffd4 = diff4 - diff3;

		if ((diffd2 < 0) && (diffd4 < 0)) {

			//
			if ((diffd2 > diffd3) && (diffd4 > diffd3)) {
				return true;
			}

		}
		return false;
	}

	private static double ave(int[] bin, int i, int j) {

		double ave = 0;
		double sum = 0;
		int cnt = 0;
		if (i < 0) {
			i = 0;
		}
		if (j > bin.length - 1) {
			j = bin.length - 1;
		}
		for (int n = i; n <= j; n++) {

			sum = sum + (bin[n]);
			cnt++;

		}
		ave = sum / (double) cnt;
		return ave;
	}

	private static int getIdx(double d) {

		double dd = d * 1000;
		int idx = (int) Math.round(dd);
		if (idx < 0)
			return 0;
		if (idx >= 4000)
			return 3999;
		return idx;
	}

}
