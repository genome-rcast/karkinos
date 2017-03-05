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
import java.util.Collection;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class CheckPossibleHDAmp {

	public static void check(DataSet dataset, PeaksInfo pi, int ploidy,
			Peak peak) throws IOException {

		if (peak == null)
			return;
		double stepsize = dataset.getBaselineLOHEstimate();

		List<CopyNumberInterval> clist = dataset.getCniVaridateList();

		float lowmean = 0;
		float highmean = 0;
		float lowestcn = 100;
		float highestcn = 0;

		List<Peak> plist = pi.getPeaklist();

		for (Peak eachpeak : plist) {

			//
			float u = (float) eachpeak.getU();

			if ((lowmean == 0) || (u < lowmean)) {
				lowmean = u;
				lowestcn = eachpeak.getCopynum();
			}
			if ((highmean == 0) || (u > highmean)) {
				highmean = u;
				highestcn = eachpeak.getCopynum();
			}

		}
		try {
			double stepsize1 = (highmean - peak.getU()) / (highestcn - ploidy);
			double stepsize2 = (peak.getU() - lowmean) / (ploidy - lowestcn);
			stepsize = Math.max(stepsize1, stepsize2);
		} catch (Exception ex) {

		}

		List<List<WaveletIF>> cap = dataset.getCapInterval();

		List<CopyNumberInterval> hdCandidates = hdCandidate(clist, lowestcn);
		List<CopyNumberInterval> ampCandidates = ampCandidate(clist, highestcn,
				cap);

		// HD
		findhdborder(hdCandidates, stepsize, cap, lowmean,lowestcn);

		// Amp
		for (CopyNumberInterval ampcand : ampCandidates) {

			// if(ampcand.getChr().equals("chr11")){
			// System.out.println("here");
			// }
			List<CopyNumberInterval> listout = checkAmp(ampcand, stepsize, cap,
					ploidy, peak.getU(), highmean, highestcn);
			if (listout != null && listout.size() > 0) {
				replace(clist, ampcand, listout);
			}
		}

	}

	private static void replace(List<CopyNumberInterval> clist,
			CopyNumberInterval ampcand, List<CopyNumberInterval> listout) {

		for (CopyNumberInterval cni : clist) {

			if (cni.equals(ampcand)) {

				clist.remove(cni);
				for(CopyNumberInterval cni2:merge(ampcand, listout)){
					if(!contain(cni2,clist)){
						clist.add(cni2);
					}
				}				
				break;

			}

		}

	}

	private static boolean contain(CopyNumberInterval cni2,
			List<CopyNumberInterval> clist) {
		
		for(CopyNumberInterval cni:clist){
			if(cni.getChr().equals(cni2.getChr())){
				if(cni.getStart() == cni2.getStart()){
					return true;
				}				
			}
		}		
		return false;
	}

	private static Collection<? extends CopyNumberInterval> merge(
			CopyNumberInterval ampcand, List<CopyNumberInterval> listout) {

		return listout;
	}

	private static List<CopyNumberInterval> checkAmp(
			CopyNumberInterval ampCand, double stepsize,
			List<List<WaveletIF>> cap, int ploidy, double u, float highmean,
			float highestcn) {

		TwoStateHMM thmm = new TwoStateHMM();

		List<WaveletIF> neighbors = getNB(ampCand, cap);
		// find steep points
		List<CopyNumberInterval> listout = thmm.checkAmp(neighbors, ampCand,
				(float) stepsize, ploidy, u, highmean, highestcn);
		//

		return listout;

	}

	private static void findhdborder(List<CopyNumberInterval> hdCandidates,
			double stepsize, List<List<WaveletIF>> cap2, float lowmean, float lowestcn)
			throws IOException {

		TwoStateHMM thmm = new TwoStateHMM();
		for (CopyNumberInterval cni : hdCandidates) {

			List<WaveletIF> neighbors = getNB(cni, cap2);
			// find steep points
			thmm.calcHD(neighbors, cni, (float) stepsize, lowmean,lowestcn);

		}

	}

	private static List<WaveletIF> getNB(CopyNumberInterval cni,
			List<List<WaveletIF>> cap2) {

		List<WaveletIF> listo = new ArrayList<WaveletIF>();

		for (List<WaveletIF> list : cap2) {

			for (WaveletIF wi : list) {
				CapInterval ci = ((CapInterval) wi);
				if (!ci.getChr().equals(cni.getChr())) {
					break;
				}

				// take wi of of 5m bp
				if (in1m(cni, ci)) {
					listo.add(wi);
				}

			}
		}

		return listo;
	}

	private static boolean in1m(CopyNumberInterval cni, CapInterval ci) {

		int unit = 1000000;
		if (cni.getStart() < ci.getEnd() + unit
				&& ci.getStart() - unit < cni.getEnd()) {
			return true;
		}
		return false;
	}

	private static List<CopyNumberInterval> ampCandidate(
			List<CopyNumberInterval> clist, float highestcn,
			List<List<WaveletIF>> cap) {
		List<CopyNumberInterval> list = new ArrayList<CopyNumberInterval>();
		for (CopyNumberInterval cni : clist) {
			boolean hit = false;
//			if (cni.getCopynumber() >= highestcn - 2) {
				if (!cni.isAllelic()) {
					list.add(cni);
					hit = true;
				}
//			}
//			if (hit == false) {
//
//				int count = 0;
//				for (List<WaveletIF> wlist : cap) {
//
//					for (WaveletIF wi : wlist) {
//
//						//
//						CapInterval ci = ((CapInterval) wi);
//						if (!ci.getChr().equals(cni.getChr())) {
//							break;
//						}
//						if (ci.getDenioseValue() >= 3.0) {
//							count++;
//						}
//
//					}
//
//				}
//				if (count >= 2) {
//					list.add(cni);
//				}
//			}
		}
		return list;

	}

	private static List<CopyNumberInterval> hdCandidate(
			List<CopyNumberInterval> clist, float lowestcn) {

		List<CopyNumberInterval> list = new ArrayList<CopyNumberInterval>();
		if (lowestcn <= 2) {
			for (CopyNumberInterval cni : clist) {

				if (cni.getCopynumber() == lowestcn) {
					list.add(cni);
				}
			}
		}
		return list;

	}

}
