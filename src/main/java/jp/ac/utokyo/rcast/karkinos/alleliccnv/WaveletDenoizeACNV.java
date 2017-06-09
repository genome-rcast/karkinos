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
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.wavelet.DistMedian;
import jwave.Transform;
import jwave.transforms.FastWaveletTransform;
import jwave.transforms.wavelets.haar.Haar1;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class WaveletDenoizeACNV {

	public static void calc(List<List<SNVHolderPlusACnv>> plist)
			throws IOException {

		int denoiselevel = getDenoiseLevel(plist) - 3;
		if (denoiselevel < 3) {
			denoiselevel = 3;
		}
		for (List<SNVHolderPlusACnv> list : plist) {
			//_calc(list, denoiselevel);
			_calcMovingAverage(list, denoiselevel);
		}
		//
		DistMedian high = new DistMedian(0.5,1.5,true);
		DistMedian low = new DistMedian(0.5,1.5,true);
		for (List<SNVHolderPlusACnv> list : plist) {
			for (SNVHolderPlusACnv sc : list) {
				double highv = sc.getHighera().getWtval();
				double lowv = sc.getLowera().getWtval();
				high.addValue(highv);
				low.addValue(lowv);
			}
		}
		System.out.println("high median =" + high.getDistributionMedian());
		System.out.println("low median =" + low.getDistributionMedian());
		double adjusth = high.getDistributionMedian() - 1;
		double adjustl = low.getDistributionMedian() - 1;
		//
		
		// adjust average
		double lowremove = 0;
		boolean norev = true;
		double mindiff  = 100;
		for (List<SNVHolderPlusACnv> list : plist) {
			for (SNVHolderPlusACnv sc : list) {

				double highv = sc.getHighera().getWtval();
				double lowv = sc.getLowera().getWtval();

				// remove value which moves in same direction
				// up to 0.2
				double highad = highv - adjusth;
				double lowad = lowv - adjustl;				
				if(lowad>highad){
					double diff = (lowad-highad);
					if(diff>lowremove){
						lowremove = diff;						
					}
					norev = false;
					if(diff<mindiff){
						diff = mindiff;
					}
				}else{
					
					double diff = (highad-lowad);
					if(diff<mindiff){
						diff = mindiff;
					}
					
				}
				
			}
		}
		
		if(mindiff==100){
			mindiff=0;
		}
		// adjust average
		for (List<SNVHolderPlusACnv> list : plist) {
			for (SNVHolderPlusACnv sc : list) {

				double highv = sc.getHighera().getWtval();
				double lowv = sc.getLowera().getWtval();

				// remove value which moves in same direction
				// up to 0.2
				double highad = highv - adjusth;
				double lowad = lowv - (adjustl + lowremove);	
				if(norev){
					lowad = (lowv - adjustl) + mindiff;
				}				
				sc.getHighera().setWtval((float) (highad));
				sc.getLowera().setWtval((float) (lowad));
			}
		}
		//
		boolean continuos = false;
		int total=0;
		int cnt =0;
		SummaryStatistics ss = new SummaryStatistics();
		for (List<SNVHolderPlusACnv> list : plist) {
			for (SNVHolderPlusACnv sc : list) {

				double highv = sc.getHighera().getWtval();
				double lowv = sc.getLowera().getWtval();
				total++;				
				if(highv<lowv){
					if(continuos){
						cnt++;
						ss.addValue(Math.abs(lowv-highv));
					}	
					
					continuos = true;
				}else{
					continuos = false;
				}

			}
		}
		//correction if low value baseline is higher than highval
		float rratio = (float)(double)((double)cnt/(double)total);
		float readjust = 0;
		if(rratio>0.2){
			readjust = (float) ss.getMean();
		}		
		//remove same directinal draft
		for (List<SNVHolderPlusACnv> list : plist) {
			for (SNVHolderPlusACnv sc : list) {

				double highv = sc.getHighera().getWtval();
				double lowv = sc.getLowera().getWtval();
				double highad = highv + readjust;
				double lowad = lowv; 
				
				boolean bothhigh = (highad > 1 && lowad > 1);
				boolean bothlow = (highad < 1 && lowad < 1);
				double denoizethre = 0.1;
				if (bothhigh || bothlow) {

					// same direction;
					double dhigh = Math.abs(1 - highad);
					double dlow = Math.abs(1 - lowad);
					double common = Math.min(dhigh, dlow);
					if (common > denoizethre) {
						common = denoizethre;
					}
					if (bothhigh) {
						highad = highad - common;
						lowad = lowad - common;
					} else {
						highad = highad + common;
						lowad = lowad + common;
					}

				}
				sc.getHighera().setWtval((float) highad);
				sc.getLowera().setWtval((float) lowad);
			}
		}
	}

	private static int getDenoiseLevel(List<List<SNVHolderPlusACnv>> plist) {

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

	private static double[] getSamplingData(
			List<List<SNVHolderPlusACnv>> plist, int samplingsize) {

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
			for (List<SNVHolderPlusACnv> list : plist) {

				if (cnt >= samplingsize)
					break;
				for (SNVHolderPlusACnv wif : list) {
					if (cnt >= samplingsize)
						break;
					testsampling[cnt] = wif.lowera.row;
					ss.addValue(wif.lowera.row);
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

	private static double lohEstimate(List<List<SNVHolderPlusACnv>> dlist) {

		float startval = 0.45f;
		float endval = 0.9f;
		double loh = lohEstimate(dlist, startval, endval);
		if (loh == startval) {
			loh = endval;
		}
		return loh;
	}

	private static double lohEstimate(List<List<SNVHolderPlusACnv>> dlist,
			float startval, float endval) {
		float maxf = 0;
		double minsd = 10;
		for (float f = startval; f <= endval; f = f + 0.01f) {

			SummaryStatistics ss = new SummaryStatistics();
			for (List<SNVHolderPlusACnv> list : dlist) {

				for (SNVHolderPlusACnv wif : list) {
					//
					double val = wif.lowera.wtval;
					if (val <= 0.9) {

						double diff = Math.abs(f - val);
						ss.addValue(diff);

					} else if (val >= 1.1) {

						double interval = Math.abs(1 - f);
						double nearlistline = getNearistLine(interval, val);
						double diff = Math.abs(nearlistline - val);
						ss.addValue(diff);

					}

				}

			}
			if (minsd > ss.getStandardDeviation()) {
				minsd = ss.getStandardDeviation();
				// System.out.println(f + "\t" + ss.getStandardDeviation());
				maxf = f;
			}
		}
		return maxf;
	}

	private static double getNearistLine(double interval, double val) {
		// up to 3n copy number gain
		double l1 = 1 + interval;
		double l2 = 1 + 2 * interval;
		//
		double diff1 = Math.abs(l1 - val);
		double diff2 = Math.abs(l2 - val);

		return diff1 < diff2 ? l1 : l2;
	}
	
	public static void _calcMovingAverage(List<SNVHolderPlusACnv> dlist, int denoiseLevel)
	throws IOException {
		
		int size = (int) Math.pow(2, denoiseLevel);
		for (int n = 0; n < dlist.size(); n++) {
			
			
			double[] ave = getMA(n,size,dlist);	
			dlist.get(n).highera.wtval = (float) ave[0];
			dlist.get(n).lowera.wtval = (float) ave[1];			
			
		}	
		
	}
	
	private static double[] getMA(int idx,int size, List<SNVHolderPlusACnv> dlist) {
		//
		double sumall0 = 0;
		double sumall1 = 0;
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
			
			SNVHolderPlusACnv  holder = dlist.get(n);
			double localweight = 1;
//			try{
//				localweight = Math.sqrt(holder.snv.getNormal().getTotal());
//			}catch(Exception ex){}	
			
			sumall0 = sumall0+(localweight*holder.highera.gcadjusted);
			sumall1 = sumall1+(localweight*holder.lowera.gcadjusted);
			weight = weight+localweight;		
		}	
		if(weight<=0)weight=1;
		double ave0 = sumall0/weight;
		double ave1 = sumall1/weight;
		return new double[]{ave0,ave1};
		
	}	

	public static void _calc(List<SNVHolderPlusACnv> dlist, int denoiseLevel)
			throws IOException {

		// extend data until data size become 2^x
		int x = 2;
		while (x < dlist.size()) {
			x = x * 2;
		}
		double[] input1 = new double[x];
		double[] input2 = new double[x];
		boolean isFirst = true;
		double setd1 = 1d;
		double setd2 = 1d;

		for (int n = 0; n < x; n++) {
			if (n < dlist.size()) {
				input1[n] = dlist.get(n).highera.gcadjusted;
				input2[n] = dlist.get(n).lowera.gcadjusted;
			} else {
				if (isFirst) {
					isFirst = false;
					setd1 = getMeanFromlast(input1, n, denoiseLevel);
					setd2 = getMeanFromlast(input2, n, denoiseLevel);
				}
				// add value to perform wavelet transform
				input1[n] = setd1;
				input2[n] = setd2;
			}
		}

		Transform t = new Transform(new FastWaveletTransform(new Haar1()));
		// System.out.println("input size="+input.length);
		double[] arrHilb1 = t.forward(input1); // 1-D AED FWT Haar forward
		double[] arrHilb2 = t.forward(input2);

		// remove small coeffictent
		int size = dlist.size() / (int) Math.pow(2, denoiseLevel);
		int nn = 0;
		for (double d : arrHilb1) {
 
			if (nn > size) {
				arrHilb1[nn] = 0;
				arrHilb2[nn] = 0;
			}
			nn++;
		}

		double[] arrReco1 = t.reverse(arrHilb1); // 1-D AED FWT Haar reverse
		double[] arrReco2 = t.reverse(arrHilb2);

		int m = 0;

		// adjustBoundary(dataList);
		for (int n = 0; n < x; n++) {
			if (n < dlist.size()) {

				dlist.get(n).highera.wtval = (float) arrReco1[n];
				dlist.get(n).lowera.wtval = (float) arrReco2[n];

			}
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

}
