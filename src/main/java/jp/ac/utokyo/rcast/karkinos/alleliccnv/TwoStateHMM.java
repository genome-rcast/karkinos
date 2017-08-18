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
import java.util.List;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;

public class TwoStateHMM {

	public static void calcHD(List<WaveletIF> neighbors,
			CopyNumberInterval cni, float stepsize, float lowmean,
			float lowestcn) throws IOException {

		Hmm<ObservationReal> hmm = getLOWHMM(neighbors, cni, stepsize, lowmean,
				lowestcn);
		List<ObservationReal> neighborsl = getList(neighbors);
		int[] hmmary = null;
		try{
		 
			hmmary = hmm.mostLikelyStateSequence(neighborsl);
		
		} catch (Exception ex) {
			return;
		}		
			
		int n = 0;

		//
		boolean startunset = true;
		int start = 0;
		int end = 0;
		for (int state : hmmary) {

			WaveletIF wi = neighbors.get(n);
			CapInterval ci = (CapInterval) wi;
			if (isStart(cni, ci) && startunset) {
				start = adjustStart(n, neighbors, hmmary);
				cni.setStart(start);
				startunset = false;
			} else if (isEnd(cni, ci)) {
				end = adjustEnd(n, neighbors, hmmary);
				cni.setEnd(end);
				break;
			}
			n++;
		}

		//
		for (WaveletIF wi : neighbors) {

			CapInterval ci = (CapInterval) wi;
			if ((ci.getStart() <= end) && (start <= ci.getEnd())) {
				ci.setCN(0);
			}
		}

	}

	private static int adjustEnd(int n, List<WaveletIF> neighbors, int[] hmmary) {

		WaveletIF wi = neighbors.get(n);
		CapInterval ci = (CapInterval) wi;

		int state = hmmary[n];
		int end = ci.getEnd();

		if (state == 0) {
			while (n < hmmary.length) {
				//
				n++;
				if (outofrange(neighbors, n)) {
					break;
				}
				wi = neighbors.get(n);
				ci = (CapInterval) wi;
				state = hmmary[n];
				if (state == 0) {
					end = ci.getEnd();
				} else {
					break;
				}

			}
		} else {
			while (n >= 0) {
				//
				n--;
				if (outofrange(neighbors, n)) {
					break;
				}
				wi = neighbors.get(n);
				ci = (CapInterval) wi;
				state = hmmary[n];
				if (state == 0) {
					end = ci.getEnd();
				} else {
					break;
				}

			}

		}
		return end;
	}

	private static int adjustStart(int n, List<WaveletIF> neighbors,
			int[] hmmary) {

		WaveletIF wi = neighbors.get(n);
		CapInterval ci = (CapInterval) wi;

		int state = hmmary[n];
		int start = ci.getStart();
		if (state == 0) {
			while (n >= 0) {
				//
				n--;
				if (outofrange(neighbors, n)) {
					break;
				}
				wi = neighbors.get(n);
				ci = (CapInterval) wi;
				state = hmmary[n];
				if (state == 0) {
					start = ci.getStart();
				} else {
					break;
				}

			}
		} else {
			while (n < hmmary.length) {
				//
				n++;
				if (outofrange(neighbors, n)) {
					break;
				}
				wi = neighbors.get(n);
				ci = (CapInterval) wi;
				state = hmmary[n];
				if (state == 0) {
					start = ci.getStart();
				} else {
					break;
				}

			}
		}
		return start;
	}

	private static boolean outofrange(List<WaveletIF> list, int n) {

		if (n < 0)
			return true;
		if (n >= list.size())
			return true;

		return false;
	}

	private static boolean isEnd(CopyNumberInterval cni, CapInterval ci) {

		int end = cni.getEnd();
		boolean bool1 = ci.getStart() <= end;
		boolean bool2 = ci.getEnd() >= end;
		return bool1 && bool2;
	}

	private static boolean isStart(CopyNumberInterval cni, CapInterval ci) {

		int start = cni.getStart();
		boolean bool1 = ci.getStart() <= start;
		boolean bool2 = ci.getEnd() >= start;
		return bool1 && bool2;
	}

	private static List<ObservationReal> getList(List<WaveletIF> neighbors) {
		List<ObservationReal> olist = new ArrayList<ObservationReal>();
		for (WaveletIF wi : neighbors) {

			double d = ((CapInterval) wi).getDenioseValue();
			// double d2 = ((CapInterval) wi).getOriginalValue();
			// double max = Math.max(d, d2);
			// ObservationReal oi = new ObservationReal(max);
			ObservationReal oi = new ObservationReal(d);
			olist.add(oi);
		}
		return olist;
	}

	private static Hmm<ObservationReal> getLOWHMM(List<WaveletIF> neighbors,
			CopyNumberInterval cni, double stepsize, float lowmean,
			float lowestcn) {

		OpdfGaussianFactory factory = new OpdfGaussianFactory();
		Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(2, factory);

		SummaryStatistics ss = getCoreVariance(neighbors, cni);
		double mean = ss.getMean();
		double variance = ss.getVariance();
		if (ss.getN() <= 1 || variance <= 0) {
			//
			if (ss.getN() == 0)
				mean = 1;
			variance = 0.5;
		}
		float minus = (float) (stepsize * 1.2);
		if (lowestcn == 2.0f) {
			minus = (float) (stepsize * 2.0);
		}
		hmm.setOpdf(0, new OpdfGaussian(lowmean - minus, variance));
		hmm.setOpdf(1, new OpdfGaussian(lowmean, variance));

		hmm.setPi(0, 0.5);
		hmm.setPi(1, 0.5);

		hmm.setAij(0, 1, 0.1);
		hmm.setAij(0, 0, 0.9);

		hmm.setAij(1, 0, 0.1);
		hmm.setAij(1, 1, 0.9);

		return hmm;
	}

	private static SummaryStatistics getCoreVariance(List<WaveletIF> neighbors,
			CopyNumberInterval cni) {
		SummaryStatistics ss = new SummaryStatistics();
		for (WaveletIF wi : neighbors) {

			CapInterval ci = (CapInterval) wi;
			int s = ci.getStart();
			int e = ci.getEnd();
			if (s < cni.getEnd() && cni.getStart() < e) {
				if (ci.getDenioseValue() > 3) {
					ss.addValue(ci.getDenioseValue());
				}
			}
		}
		if (ss.getN() == 0) {
			for (WaveletIF wi : neighbors) {
				CapInterval ci = (CapInterval) wi;
				ss.addValue(ci.getDenioseValue());
			}
		}
		return ss;
	}

	public List<CopyNumberInterval> checkAmp(List<WaveletIF> neighbors,
			CopyNumberInterval cni, float stepsize, int ploidy, double u,
			float highmean, float highestcn) {

		//
		SummaryStatistics ss = getCoreVariance(neighbors, cni);
		//
		float copygain = (float) ((ss.getMean() - u) / stepsize) + ploidy;
//		System.out.println(cni.getChr() +" copygain=" + copygain);
//		if (copygain <= 2.0) {
//			return null;
//		}
		// if (copygain <= 3.5) {
		// return null;
		// }

		Hmm<ObservationReal> hmm = getHighHMM(neighbors, cni, stepsize, ss,
				highmean, highestcn);

		List<ObservationReal> neighborsl = getList(neighbors);
		int[] hmmary = null;
		try {
			hmmary = hmm.mostLikelyStateSequence(neighborsl);
		} catch (Exception ex) {
			return null;
		}

		List<CopyNumberInterval> list = new ArrayList<CopyNumberInterval>();

		//
		int n = 0;
		CopyNumberInterval cnin = null;
		CapInterval cib4 = null;
		SummaryStatistics ss2 = new SummaryStatistics();
		int cnt = 0;
		for (int state : hmmary) {

			WaveletIF wi = neighbors.get(n);
			CapInterval ci = (CapInterval) wi;
			if (state == 1) {
				cnt++;
				ss2.addValue(ci.getValue());
			}
			n++;
		}
		if (cnt > 0) {

		} else {
			return null;
		}
		copygain = (float) ((ss2.getMean() - u) / stepsize) + ploidy;
		if (copygain <= (ploidy*2)) {
			return null;
		}
		n = 0;
		int stateb4 = -1;
		for (int state : hmmary) {

			WaveletIF wi = neighbors.get(n);
			CapInterval ci = (CapInterval) wi;
			if(state==0){
				if(ci.getDenioseValue() > 3.5){
					state = 1;
				}
			}			
			
			boolean statechange = (stateb4 != state);

			if (statechange) {

				if (state == 0) {

					if (cnin != null) {
						list.add(cnin);
						cnin = null;
					}
					cnin = new CopyNumberInterval();
					cnin.setChr(ci.getChr());
					cnin.setStart(ci.getStart());
					cnin.setCopynumber(cni.getCopynumber());

				} else {

					if (cnin != null) {
						list.add(cnin);
						cnin = null;
					}
					cnin = new CopyNumberInterval();
					cnin.setHdelation(true);
					cnin.setChr(ci.getChr());
					cnin.setStart(ci.getStart());
					cnin.setCopynumber(copygain);
				}

			} else {
				if (cnin != null) {
					cnin.setEnd(ci.getEnd());

				}
			}

			n++;
			stateb4 = state;
		}

		if (cnin != null) {
			list.add(cnin);
		}

		return list;
	}

	private Hmm<ObservationReal> getHighHMM(List<WaveletIF> neighbors,
			CopyNumberInterval cni, float stepsize, SummaryStatistics ss,
			float highmean, float highestcn) {

		OpdfGaussianFactory factory = new OpdfGaussianFactory();
		Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(2, factory);

		//
		double mean = ss.getMean();
		double variance = ss.getVariance();
		if (ss.getN() <= 1 || variance <= 0) {
			//
			if (ss.getN() == 0)
				mean = 1;
			variance = 0.5;
		}
		float highmean2 = highmean + (2 * stepsize);
		if (highmean2 > 3.0) {
			highmean2 = 3.0f;
		}
		hmm.setOpdf(1, new OpdfGaussian(highmean2, variance));
		hmm.setOpdf(0, new OpdfGaussian(highmean, variance));

		hmm.setPi(0, 0.5);
		hmm.setPi(1, 0.5);

		hmm.setAij(0, 1, 0.1);
		hmm.setAij(0, 0, 0.9);

		hmm.setAij(1, 0, 0.1);
		hmm.setAij(1, 1, 0.9);
		return hmm;

	}

}
