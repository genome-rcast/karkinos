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
import java.util.Collections;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jwave.Transform;
import jwave.transforms.FastWaveletTransform;
import jwave.transforms.wavelets.haar.Haar1;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class WaveletDenoize {

	public static double calc(DataSet dataset) throws IOException {

		List<List<WaveletIF>> plist = dataset.getCapInterval();
		int denoiselevel = getDenoiseLevel(plist);
		for (List<WaveletIF> list : plist) {
			//_calc(list, denoiselevel);
			//
			_calcMovingAverage(list, denoiselevel);
		}
		double baselineLOHEstimate = lohEstimate(plist);
		dataset.setBaselineLOHEstimate(baselineLOHEstimate);
		for (List<WaveletIF> list : plist) {
			setcn(list, baselineLOHEstimate);
		}
		return baselineLOHEstimate;
		// System.out.println(firstEstimateTumorRate);

	}

	private static int getDenoiseLevel(List<List<WaveletIF>> plist) {

		int samplingsize = 4096;
		int denoizemax = KarkinosProp.maxdenoiseLevel;
		double[] testsampling = getSamplingData(plist, samplingsize);
		if (testsampling == null) {
			return denoizemax;
		}

		int n = 1;
		for (; n < denoizemax; n++) {

			SummaryStatistics ss = getSD(testsampling);
			double sd = ss.getStandardDeviation();
			System.out.println("mean=" + ss.getGeometricMean() + "sd=" + sd);
			if (sd < KarkinosProp.denozeToSD)
				break;
			testsampling = downsampling(testsampling);

		}
		System.out.println("denoise level = " + (n));
		return n;

	}

	private static double[] getSamplingData(List<List<WaveletIF>> plist,
			int samplingsize) {

		double minmeandiff = 0;
		double[] testsamplingr = null;
		// take sampling data where no CNV occurring
		// try 10 times from chr1 start, then chr2 start so on
		// take minimum difference sampleset from n=1
		int m = 0;
		for (int chridx = 0; chridx < 5; chridx++) {

			SummaryStatistics ss = new SummaryStatistics();
			double[] testsampling = new double[samplingsize];
			int cnt = 0;
			if (m == chridx)
				continue;
			for (List<WaveletIF> list : plist) {

				if (cnt >= samplingsize)
					break;
				for (WaveletIF wif : list) {
					if (cnt >= samplingsize)
						break;
					testsampling[cnt] = wif.getValue();
					ss.addValue(wif.getValue());
					cnt++;
				}

			}
			if (minmeandiff == 0) {
				minmeandiff = (1 - ss.getMean());
				testsamplingr = testsampling;
			}
			if (minmeandiff > (1 - ss.getMean())) {
				minmeandiff = (1 - ss.getMean());
				testsamplingr = testsampling;
			}
			m++;
		}
		return testsamplingr;
	}

	private static double[] downsampling(double[] testsampling) {

		double[] newa = new double[testsampling.length / 2];
		for (int n = 0; n + 1 < testsampling.length; n = n + 2) {

			double d0 = testsampling[n];
			double d1 = testsampling[n + 1];
			double mean = (d0 + d1) / 2;
			newa[n / 2] = mean;
		}
		return newa;
	}

	private static SummaryStatistics getSD(double[] testsampling) {
		SummaryStatistics ss = new SummaryStatistics();
		for (double d : testsampling) {
			ss.addValue(d);
		}
		return ss;
	}

	private static void setcn(List<WaveletIF> list, double baselineLOHEstimate) {

		for (WaveletIF wif : list) {

			double denoise = wif.getDenioseValue();
			double nval = step(denoise, baselineLOHEstimate);
			wif.setCN(nval);

		}
		adjustBoundary(list);

	}

	private static double lohEstimate(List<List<WaveletIF>> dlist) {
		
		//interval0
		boolean detected0 = true;
		float startval0 = 0.45f;
		float endval0 = 0.75f;
		double loh0 = lohEstimate(dlist, startval0, endval0,false);
		if (loh0 == startval0||loh0 == endval0) {
			loh0 = endval0;
			detected0 = false;
		}
		
		//interval1	
		
		boolean detected = true;
		float startval = 0.40f;
		float endval = 0.9f;
		double loh = lohEstimate(dlist, startval, endval,true);		
		if (loh == startval||loh0 == endval) {
			loh = endval;
			detected = false;
		}
		
		float margin = 0.05f;
		
		if(detected0&&detected){
			
			boolean morelwpeak = moreThanParcentLowArea(0.02f,loh-margin, dlist);
			if(morelwpeak){
				return Math.min(loh0,loh);
			}else{
				return Math.max(loh0,loh);
			}
			
		}else if(detected){
			return loh;
		}else if(detected0){
			return loh0;
		}		
		return loh;
	}

	private static double lohEstimate(List<List<WaveletIF>> dlist,
			float startval, float endval,boolean include3n) {
		float maxf = 0;
		float maxf2 = 0;
		double minsd = 10;
		double minsd2 = 10;

		List<double[]> l = new ArrayList<double[]>();
		List<double[]> l_low = new ArrayList<double[]>();
		
		for (float f = startval; f <= endval; f = f + 0.01f) {

			SummaryStatistics ss = new SummaryStatistics();
			SummaryStatistics sslow = new SummaryStatistics();
						
			for (List<WaveletIF> list : dlist) {

				for (WaveletIF wif : list) {
					//
					double val = wif.getDenioseValue();
					if (val <= endval) {

						double diff = Math.abs(f - val);
						ss.addValue(diff);
						sslow.addValue(diff);

					} else if (val >= 1+Math.abs(1-endval)) {

						double interval = Math.abs(1 - f);
						double nearlistline = getNearistLine(interval, val,include3n);
						//double nearlistline = 1 + interval;
						double diff = Math.abs(nearlistline - val);
						if(diff<interval){
							ss.addValue(diff);
						}
						
					}

				}

			}
			//System.out.println(f + "\t" + ss.getStandardDeviation());
			if (minsd > ss.getStandardDeviation()) {

				minsd = ss.getStandardDeviation();
				System.out.println(f + "\t" + ss.getStandardDeviation());
				maxf = f;
			}
			l.add(new double[] { f, ss.getStandardDeviation() });
			l_low.add(new double[] { f, sslow.getStandardDeviation() });
			

		}
		float lowpeak = getLowerPeak(l,l_low,dlist);
		if (lowpeak > 0) {
			maxf = lowpeak;
		}
		return maxf;
	}

	private static float getLowerPeak(List<double[]> l,List<double[]> l_low,
			List<List<WaveletIF>> dlist) {

		List<double[]> peaks = new ArrayList<double[]>();
		for (int n = 0; n + 2 < l.size(); n++) {

			double f0 = l.get(n)[0];
			double sd0 = l.get(n)[1];

			double f1 = l.get(n + 1)[0];
			double sd1 = l.get(n + 1)[1];

			double f2 = l.get(n + 2)[0];
			double sd2 = l.get(n + 2)[1];
						
			if ((sd0 > sd1) && (sd2 > sd1)) {
				
				peaks.add(l.get(n + 1));
			}

		}
			

		if (peaks.size() <= 1) {
			return -1;
		} else {

			//
			sort(peaks);
			double[] firstPeak = peaks.get(0);
			double[] secondPeak = peaks.get(1);
			if (firstPeak[0] < secondPeak[0]) {
				return (float) firstPeak[0];
			} else {
				float parcent = 0.03f;
				float margin = 0.05f;
				if (firstPeak[0] - margin> secondPeak[0]) {
					if (moreThanParcentLowArea(parcent, firstPeak[0]-margin, dlist)) {

						return (float) secondPeak[0];

					}
				}

			}
			return (float) firstPeak[0];
		}

	}

	private static boolean moreThanParcentLowArea(float f, double firstPeak,
			List<List<WaveletIF>> dlist) {
		try {
			long total = 0;
			long lowcounts = 0;
			for (List<WaveletIF> list : dlist) {

				for (WaveletIF wif : list) {
					//
					double val = wif.getDenioseValue();
					total++;
					if (val < firstPeak) {
						lowcounts++;
					}
//					if (val > 1+(1-firstPeak)) {
//						lowcounts++;
//					}
					
				}
			}
			double ratio = (double) lowcounts / (double) total;
			return ratio > f;
		} catch (Exception ex) {

		}
		return false;
	}

	private static void sort(List<double[]> peaks) {

		Collections.sort(peaks, new MyCompPeak());

	}

	private static double getNearistLine(double interval, double val, boolean include3n) {
		// up to 3n copy number gain
		double l1 = 1 + interval;
		if(include3n=false){
			return l1;
		}		
		double l2 = 1 + 2 * interval;
		//
		double diff1 = Math.abs(l1 - val);
		double diff2 = Math.abs(l2 - val);

		return diff1 < diff2 ? l1 : l2;
	}

	public static void _calcMovingAverage(List<WaveletIF> dlist, int denoiseLevel){
		
		//
		int size = (int) Math.pow(2, denoiseLevel);
		for (int n = 0; n < dlist.size(); n++) {
			
			
			WaveletIF wif = dlist.get(n);			
			double ave = getMA(n,size,dlist);	
			if(Double.isNaN(ave)){
				ave = wif.getValue();
			}			
			wif.setDenioseValue(ave);			
			
		}		
	
	}
	
	private static double getMA(int idx,int size, List<WaveletIF> dlist) {
		//
		double sumall = 0;
		double sumbase = 0;
		double weight = 0;
		int half = size/2;
		int start = idx-half;
		int addend = 0;
		if(start<0){
			addend = Math.abs(start);
			start =0;
		}
		int end = idx + half+addend;
		if(end>=dlist.size()){
			
			int minusstart = end-dlist.size();
			start = start - minusstart;
			if(start<0){
				start = 0;
			}
			end = dlist.size();
		}
		
		for(int n=start;n<end;n++){
			
			WaveletIF wi = dlist.get(n);
			CapInterval ci = (CapInterval)wi;
			//sumbase = sumbase + ci.getLength();
			double localweight = 1;
			try{				
				localweight = Math.sqrt(ci.getCNVInfo().getNormalcnt());
			}catch(Exception ex){}			
			weight  = weight + localweight;
			
			double val = wi.getValue();
			//sumall = sumall+(val*ci.getLength());
			sumall = sumall+(val*localweight);
		}		
		double ave = sumall/weight;
		return ave;
		
	}	

	public static void _calc(List<WaveletIF> dlist, int denoiseLevel)
			throws IOException {

		// extend data until data size become 2^x
		int x = 2;
		while (x < dlist.size()) {
			x = x * 2;
		}
		double[] input = new double[x];
		boolean isFirst = true;
		double setd = 1d;
		for (int n = 0; n < x; n++) {
			if (n < dlist.size()) {
				input[n] = dlist.get(n).getValue();
			} else {
				if (isFirst) {
					isFirst = false;
					setd = getMeanFromlast(input, n, denoiseLevel);
				}
				// add value to perform wavelet transform
				input[n] = setd;
			}
		}

		Transform t = new Transform(new FastWaveletTransform(new Haar1()));
		// System.out.println("input size="+input.length);
		double[] arrHilb = t.forward(input); // 1-D AED FWT Haar forward

		// remove small coeffictent
		int size = dlist.size() / (int) Math.pow(2, denoiseLevel);
		int nn = 0;
		for (double d : arrHilb) {

			if (nn > size) {
				arrHilb[nn] = 0;
			}
			nn++;
		}

		double[] arrReco = t.reverse(arrHilb); // 1-D AED FWT Haar reverse

		int m = 0;
		List<Double[]> dataList = new ArrayList<Double[]>();
		for (double denoiseval : arrReco) {

			// System.out.println(d);
			WaveletIF wif = dlist.get(m);
			// double copynumber = step(denoiseval);
			Double[] da = new Double[3];
			da[0] = wif.getValue();
			da[1] = denoiseval;
			// da[2]=copynumber;
			dataList.add(da);
			m++;
			if (m > dlist.size() - 1) {
				break;
			}

		}
		// adjustBoundary(dataList);
		m = 0;
		for (Double[] da : dataList) {
			WaveletIF wif = dlist.get(m);
			wif.setDenioseValue(da[1]);
			// wif.setCN(da[2]);
			m++;
		}

	}

	private static double getMeanFromlast(double[] input, int n,
			int denoiselevel) {

		int cnt = 0;
		int samplesize = (int) Math.pow(2, denoiselevel);
		SummaryStatistics ss = new SummaryStatistics();
		for (int m = n; m > 0 && cnt < samplesize; m--) {
			cnt++;
			ss.addValue(input[m]);
		}
		return ss.getMean();
	}

	private static void adjustBoundary(List<WaveletIF> list) {

		int m = 0;
		double cnb4 = 0;
		for (WaveletIF wif : list) {

			double copynumber = wif.getCN();
			if (cnb4 == 0) {
				cnb4 = copynumber;
			}
			//
			if (cnb4 != copynumber) {
				boundCheck(list, m);
			}
			cnb4 = copynumber;
			m++;
		}

	}

	private static void boundCheck(List<WaveletIF> list, int n) {

		int len = list.size();
		int m = n;
		while (m > 0) {

			WaveletIF da1 = list.get(m);
			double copynumber1 = da1.getCN();
			WaveletIF da0 = list.get(m - 1);
			double row0 = da0.getValue();
			double copynumber0 = da0.getCN();
			if (Math.abs(copynumber1 - row0) < Math.abs(copynumber0 - row0)) {
				da0.setCN(copynumber1);
			} else {
				break;
			}
			m--;
		}
		int m2 = n - 1;
		while (m2 + 1 < len) {

			WaveletIF da1 = list.get(m2);
			double copynumber1 = da1.getCN();
			WaveletIF da2 = list.get(m2 + 1);
			double row2 = da2.getValue();
			double copynumber2 = da2.getCN();
			if (Math.abs(copynumber1 - row2) < Math.abs(copynumber2 - row2)) {
				da2.setCN(copynumber1);
			} else {
				break;
			}
			m2++;
		}

	}

	private static double step(double d, double lohbase) {
		// simple stem function
		float stepd = (float) lohbase;
		float stepsize = (float) (1 - lohbase);
		double mindiff = 0;
		int out = 0;
		int cnt = 0;
		for (; stepd <= 5; stepd = stepd + stepsize) {

			double diff = Math.abs(stepd - d);
			if (mindiff == 0) {
				mindiff = diff;
				out = cnt;
			}
			if (diff < mindiff) {
				mindiff = diff;
				out = cnt;
			}
			cnt++;
		}
		return (out + 1) * 0.5;
	}

}
