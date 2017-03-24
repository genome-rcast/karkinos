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

	public static double vinit = 0.0003;

	public static PeaksInfo calc(DataSet dataset, double lohestimate)
			throws IOException {

		PeaksInfo pi = new PeaksInfo();
		// take over all bin distribution
		double dmin = 0;
		double dmax = 4;
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
				int idx = getIdx(val);
				bininit[idx] = bininit[idx] + 1;
				if (Math.abs(1 - val) < 0.25) {
					ss.addValue(val);
				}

			}

		}
		double aved = total / 4000;
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
		double prevmin = 0;
		double prevminy = 0;
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
					prevmin = 0;
					prevminy = 0;
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
		int cntflg = 0;
		List<EMval> list = new ArrayList<EMval>();
		for (List<WaveletIF> dlist : plist) {

			for (int n = 0; n < dlist.size(); n++) {

				WaveletIF wif = dlist.get(n);
				EMval emval = new EMval();
				emval.setValue(wif.getValue());
				list.add(emval);
			}

		}

		// EM not working well
		// excute em
		// int count = 0;
		// while (cntflg == 0) {
		//
		// estep(list, peaklist);
		// int ret = mstep(list, peaklist);
		// if (ret < 0) {
		// break;
		// }
		// if (count > 1) {
		// break;
		// }
		// count++;
		// }
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
		int maxpeakidx = 0;
		double maxr = 0;
		double minu = 0;
		for (Peak p : peaklist) {
			double r = p.getR() / psum;
			if (r > maxr) {
				maxr = r;
				maxpeakidx = pidx;
			}
			
			if((minu==0)||(minu > p.getU())){
				minu = p.getU();
			}
			
			p.setR(r);
			p.setPeakidx(pidx);
			pidx++;
			System.out.println(p.getR() + "\t" + p.getU() + "\t" + p.getV());
		}

		// find peak distance
//		float peakdist = getUnitPeakDist(peaklist, maxpeakidx,lohestimate);
//		// add artifitial peaks 4 low 4 high
//		List<Peak> peaklistn = addArtifitialLowPeaks(minu,peaklist, maxpeakidx, peakdist);
//		pi.setPeaklist(peaklistn);
		return pi;

	}

	private static List<Peak> addArtifitialLowPeaks(double minu,List<Peak> peaklist, int maxpeakidx,
			float peakdist) {
		
		Peak maxpeak = peaklist.get(maxpeakidx);
		double maxu = maxpeak.getU();
		double sd = maxpeak.getSD()*2;
		
		List<Peak> peaklistnew 
			= new ArrayList<Peak>();
		
		for (int n = -4; n <= -1; n++) {
			
			double u = maxu+(peakdist*n);
			if(u>minu)break;
			if(notContain(peaklist,u,sd)&&(u>0)){
				
				
				Peak p = new Peak(u, 100);
				p.setU(u);
				// p.setV(ss.getVariance());
				p.setV(maxpeak.getV());
				p.setDefault(true);
				p.setR(maxpeak.getR() * 0.05);
				p.setArtifitial(true);
				p.setDefault(true);
				peaklistnew.add(p);
			}
			
		}
		peaklistnew.addAll(peaklist);
		return peaklistnew;

	}

	private static boolean notContain(List<Peak> peaklist, double u,double sd) {
	
		if(u<0)return false;
		for(Peak p:peaklist){
			
			if(Math.abs(p.getU()-u)<sd){
				return false;
			}
			
		}		
		return true;
	}

	private static float getUnitPeakDist(List<Peak> peaklist, int maxpeakidx, double lohestimate) {

		Peak maxpeak = peaklist.get(maxpeakidx);
		double sd = maxpeak.getSD();
		Map<Float, Float> dists = new TreeMap<Float, Float>();
		for (int n = 0; n < peaklist.size(); n++) {

			if (n == maxpeakidx)
				continue;
			Peak peak = peaklist.get(n);
			float dist = (float) Math.abs(peak.getU() - maxpeak.getU());
			dists.put(dist, (float) peak.getR());

		}
		List<PeakDistMagnitudeBean> list = new ArrayList<PeakDistMagnitudeBean>();
		//
		Iterator<Float> ite = dists.keySet().iterator();
		while (ite.hasNext()) {

			//
			float pd = ite.next();
			float mg = dists.get(pd);
			reg(pd, mg, sd, list);

		}

		int maxreg = 0;
		PeakDistMagnitudeBean maxbean = null;
		for (PeakDistMagnitudeBean bean : list) {
			if (bean.map.size() > maxreg) {
				maxreg = bean.map.size();
				maxbean = bean;
			}

		}
		if (maxbean == null) {
			return (float) lohestimate;
		}
		return maxbean.getDist();

	}

	private static void reg(float pd, float mg, double sd,
			List<PeakDistMagnitudeBean> list) {

		PeakDistMagnitudeBean bean = null;
		for (PeakDistMagnitudeBean bb : list) {

			//
			for (int n = 1; n <= 6; n++) {

				boolean include = include(bb.getDist(), pd / n, sd);
				if (include) {
					bean = bb;
					break;
				}

			}

		}
		if (bean == null) {
			bean = new PeakDistMagnitudeBean();
			list.add(bean);
		}
		bean.reg(pd, mg);

	}

	private static boolean include(float dist, float pd, double sd) {

		return Math.abs(dist - pd) < sd;

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

	// private static double getConvolve(double[] binma, int n) {
	//
	// double sum = 0;
	// for(int m=0;m<200;m++){
	//
	// if(m<50)continue;
	// if(m>149)continue;
	// //
	// double val = binma[(n-100)+m];
	// //convert 0 to -4.5
	// // 200 to 4.5
	// double convval = convval(m);
	// double multiple = getNormal(convval);
	// sum = sum + (val*multiple);
	//
	// }
	// return sum;
	// }

	private static double convval(int m) {

		int n = m - 100;
		return (n / 100) * 9;

	}

	private static double getNormal(double convval) {

		double nd = (1.0 / (Math.sqrt(2 * Math.PI)))
				* Math.exp(-0.5 * Math.pow(convval, 2));
		return nd;
	}

	private static List<Peak> assignCNtoPeak(List<Peak> peakslist,
			List<EMval> list, double lohestimate, double variance) {

		int idx = 0;
		int n1idx = 0;
		double n1level = 0;
		double maxratio = 0;
		// assgin n = 2
		for (Peak p : peakslist) {

			if (p.getR() > maxratio && (p.getU() > 0.75) && (p.getU() < 1.25)) {
				maxratio = p.getR();
				n1level = p.getU();
				n1idx = idx;
			}
			idx++;
		}
		Peak p2n = peakslist.get(n1idx);
		p2n.setCn(1.0f);
		double sd2n = p2n.getSD();
		//
		List<Peak> peakslistnew = new ArrayList<Peak>();
		// remove peak with too broad variance and set variance to init val
		for (Peak p : peakslist) {

			if (p.getR() > 0.03 || p.getSD() < (sd2n * 2)) {

				//
				peakslistnew.add(p);

			}

		}
		// EM with new distribution
		// int count = 0;
		// while (true) {
		//
		// estep(list, peakslistnew);
		// int ret = mstep(list, peakslistnew);
		// if (ret < 0) {
		// break;
		// }
		// if (count > 1) {
		// break;
		// }
		// count++;
		// }

		// find max ratio peak under 2n LOH
		// LOH section
		// LOH
		Peak p2n_2 = null;
		for (Peak p : peakslistnew) {

			//
			if (p.getCn() == 1.0f) {

				// normal line
				p2n_2 = p;
				p2n_2.setDefault(true);

			}

		}

		Peak p1n_2 = null;
		double r = 0;
		double stepsize = 0;
		for (Peak p : peakslistnew) {

			if (p.getCn() > 0)
				continue;
			if (p.getU() < p2n_2.getU()) {

				if (p.getR() > 0.05) {
					// lowest and hava more than 5%
					p1n_2 = p;
					stepsize = Math.abs(p.getU() - p2n_2.getU());
					break;
				} else {

					// else max exsist ratio
					if (p.getR() > 0.01
							&& withinParcent(0, 0.07, lohestimate, p.getU())) {

						if (p.getR() > r) {
							stepsize = Math.abs(p.getU() - p2n_2.getU());
							p1n_2 = p;
							r = p.getR();
						}
					}
				}

			}

		}
		//
		if (p1n_2 != null) {
			p1n_2.setCn(0.5f);
			p1n_2.setDefault(true);

		}
		if (stepsize == 0) {

			// no observed LOH
			stepsize = Math.abs(1 - lohestimate);
			// add vartual peak
			double uloh = p2n_2.getU() - stepsize;
			Peak p = new Peak(uloh, 100);
			p.setU(uloh);
			// p.setV(ss.getVariance());
			p.setV(variance);
			p.setCn(0.5f);
			p.setDefault(true);
			p.setR(p2n_2.getR() * 0.05);
			p.setArtifitial(true);
			p.setDefault(true);
			peakslistnew.add(p);

		}

		// set 2n,3n
		double twonline = p2n_2.getU();
		boolean n2set = false;
		boolean n3set = false;
		for (int n = 1; n <= 6; n++) {

			//
			double line = twonline + (n * stepsize);
			Peak passign = null;
			for (Peak p : peakslistnew) {

				if (p.getCn() > 0)
					continue;
				double val = p.getU();
				boolean fit = withinParcent(n, 0.07, line, val);
				if (fit) {

					//
					if ((passign == null) || (passign.getR() < p.getR())) {
						passign = p;
					}

				}

			}
			if (passign != null) {

				if (n == 1) {
					n2set = true;
					passign.setDefault(true);
				}
				if (n == 2) {
					n3set = true;
					passign.setDefault(true);
				}
				passign.setCn(1 + (n * 0.5f));

			}

		}
		//
		if (n2set == false) {

			double line = twonline + stepsize;
			Peak p = new Peak(line, 100);
			p.setU(line);
			// p.setV(ss.getVariance());
			p.setV(variance);
			p.setCn(1.5f);
			p.setR(p2n_2.getR() * 0.05);
			p.setArtifitial(true);
			p.setDefault(true);
			peakslistnew.add(p);

		}
		if (n3set == false) {

			double line = twonline + (2 * stepsize);
			Peak p = new Peak(line, 100);
			p.setU(line);
			// p.setV(ss.getVariance());
			p.setV(variance);
			p.setCn(2.0f);
			p.setR(p2n_2.getR() * 0.05);
			p.setArtifitial(true);
			p.setDefault(true);
			peakslistnew.add(p);

		}

		sort(peakslistnew);
		List<Peak> peakslistnew2 = new ArrayList<Peak>();
		if (variance > 0.03) {
			variance = 0.03;
		}

		// assing subpopulation peaks
		for (Peak p : peakslistnew) {

			if (p.getCn() > 0) {
				p.setV(variance);
				peakslistnew2.add(p);

			} else if (p.getCn() == 0) {

				double u = p.getU();
				if (p.getR() < 0.005) {
					continue;
				}
				if (outofmain2sd(u, peakslistnew)) {
					double cn = 0;
					if (u > twonline) {
						cn = 1 + ((Math.abs(u - twonline) / stepsize) * 0.5f);
					} else {
						cn = 1 - ((Math.abs(u - twonline) / stepsize) * 0.5f);
					}
					if (cn < 0) {
						cn = 0;
					}
					if (cn > 4)
						cn = 4;
					p.setCn((float) cn);
					p.setV(variance);
					peakslistnew2.add(p);
				}

			}

		}
		return peakslistnew2;

	}

	private static boolean outofmain2sd(double u, List<Peak> peakslist) {

		for (Peak p : peakslist) {

			if (p.getDefault()) {

				double u2 = p.getU();
				double sd2 = (p.getSD() * 2);
				if ((u < (u2 + (sd2))) && (u > (u2 - (sd2)))) {
					return false;
				}

			}
		}
		return true;
	}

	private static void sort(List<Peak> peakslistnew) {

		Collections.sort(peakslistnew, new PeakComparator());

	}

	private static boolean withinParcent(int n, double d, double val1,
			double val2) {

		double diffwithin = val1 * d;
		if (n == 0) {
			n = 1;
			diffwithin = diffwithin * 2;
		}
		if (n > 1) {
			double add = diffwithin * n * 0.5;
			diffwithin = diffwithin + add;
		}
		return Math.abs(val1 - val2) < diffwithin;

	}

	public static double getNdistP(double x, double u, double v) {

		double p = Math.exp(-0.5 * (Math.pow((x - u), 2) / v))
				/ Math.sqrt(2 * Math.PI * v);

		if (Double.isNaN(p)) {
			p = 0;
		}
		return p;
	}

	private static void estep(List<EMval> list, List<Peak> peaklist) {

		// /
		for (EMval emval : list) {

			double likesum = 0;
			for (Peak peak : peaklist) {

				//
				double r = peak.getR();
				double v = peak.getV();
				double u = peak.getU();
				double pnd = getNdistP(emval.getValue(), u, v);
				double like = r * pnd;
				likesum = likesum + like;

			}
			List<Double> peakPlist = new ArrayList<Double>();
			for (Peak peak : peaklist) {

				//
				double r = peak.getR();
				double v = peak.getV();
				double u = peak.getU();
				double pnd = getNdistP(emval.getValue(), u, v);
				double like = r * pnd;
				double p = 0;
				if (likesum > 0) {
					p = like / likesum;
				}
				peakPlist.add(p);
			}
			emval.setPforpeak(peakPlist);

		}

	}

	private static int mstep(List<EMval> list, List<Peak> peaklist) {

		for (Peak peak : peaklist) {

			peak.setSUMZ(0);
			peak.setSUMZX(0);
			peak.setSUMZX2(0);

		}

		for (EMval emval : list) {

			int idx = 0;
			for (Peak peak : peaklist) {

				double z = emval.getPforpeak().get(idx);
				double SUMZ = peak.getSUMZ() + z;
				double SUMZX = peak.getSUMZX() + (z * emval.getValue());
				peak.setSUMZ(SUMZ);
				peak.setSUMZX(SUMZX);

				idx++;
			}
		}

		// r,u
		for (Peak peak : peaklist) {

			double SUMZ = peak.getSUMZ();
			double SUMZX = peak.getSUMZX();

			double r = SUMZ / (double) list.size();
			double u = SUMZX / SUMZ;

			if (Double.isNaN(r) || Double.isNaN(u)) {
				return -1;
			}

			// renew r and u
			peak.setR(r);
			peak.setU(u);

		}

		// u
		for (EMval emval : list) {

			int idx = 0;
			for (Peak peak : peaklist) {

				double z = emval.getPforpeak().get(idx);
				double SUMZX2 = peak.getSUMZX2()
						+ (z * pow2(emval.getValue() - peak.getU()));
				peak.setSUMZX2(SUMZX2);
				idx++;
			}
		}

		for (Peak peak : peaklist) {

			double SUMZ = peak.getSUMZ();
			double SUMZX2 = peak.getSUMZX2();

			double v = SUMZX2 / SUMZ;
			if (Double.isNaN(v)) {
				return -1;
			}

			if (v < peak.getVinit() * 2) {
				peak.setV(v);
			}

		}
		return 1;
	}

	private static double pow2(double x) {
		return Math.pow(x, 2);
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
