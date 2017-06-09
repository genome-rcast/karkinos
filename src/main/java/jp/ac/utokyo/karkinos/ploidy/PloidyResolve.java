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
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.SNVHolderPlusACnv;
import jp.ac.utokyo.rcast.karkinos.distribution.DataHolderByCN;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.wavelet.ChildPeak;
import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class PloidyResolve {

	MatchMatrixBean mmb1 = new MatchMatrixBean();

	public MatchMatrixBean getMmb1() {
		return mmb1;
	}

	public void resolve(DataSet dataset, AllelicCNV alCNV, PeaksInfo pi,
			int cnvcount) {

		Map<Integer, PeakAnalysisComponent> map = new LinkedHashMap<Integer, PeakAnalysisComponent>();
		Map<SNVHolder, SNVHolderPlusACnv> hsmap = toMap(alCNV);

		//
		int t0 = 0, t1 = 0, t2 = 0;
		double d0 = 0, d1 = 0, d2 = 0;
		int pidxb4 = -1;
		int cntn = 0;
		for (SNVHolder snv : dataset.getSnvlist()) {

			//

			cntn++;

			int pidx = (int) snv.getCi().getPeakIdx();

			PeakAnalysisComponent pac = null;
			if (map.containsKey(pidx)) {
				//
				pac = map.get(pidx);

			} else {
				pac = new PeakAnalysisComponent();
				map.put(pidx, pac);

			}
			if (snv.isHetroSNP()) {

				if (pidxb4 != pidx) {
					System.out.println(cntn + "_" + snv.getPos() + "_" + pidx);
				}


				pidxb4 = pidx;
				//
				if (pidx == 0) {
					t0++;
					d0 = d0 + snv.getTumor().getRatio();
				} else if (pidx == 1) {
					t1++;
					d1 = d1 + snv.getTumor().getRatio();
				} else {
					t2++;
					d2 = d2 + snv.getTumor().getRatio();
				}

				// System.out.println(snv.getTumor().getRatio());
				pac.setSNPInfo(snv);
				if (hsmap.containsKey(snv)) {
					pac.setSNPPlus(hsmap.get(snv));
				}

			}
			int flg = snv.getFlg();
			boolean check = (flg == PileUP.SomaticMutation);
			if (check) {

				// boolean pass = snv.getFilterResult().isPassFilter();
				// boolean enoughdepth = snv.getTumor().getTotalcnt() > 50;
				// if(pass){
				// pac.setSMutaion(snv);
				// }

			}
			//System.out.println(cntn + "_" + snv.getPos() + "_" + pidx);
			
		}
		
		System.out.println(t0 + " " + t1 + " " + t2);
		System.out.println(d0 / t0 + " " + d1 / t1 + " " + d2 / t2);

		//
		Iterator<Integer> peakidxs = map.keySet().iterator();
		int maxMagnatudeEvenIdx = -1;
		double maxmagnitude = 0;
		List<PeakTmpInfo> evenpeaks = new ArrayList<PeakTmpInfo>();
		while (peakidxs.hasNext()) {

			int peakidx = peakidxs.next();
			PeakAnalysisComponent pac = map.get(peakidx);
			Peak p = pi.getPeaklist().get(peakidx);
			pac.analyse();

			if (pac.evenpeak || pac.complexpeak) {

				//
				double magnitude = p.getR();
				if (pac.complexpeak) {
					magnitude = p.getR() * pac.getEvenR();

				}
				if (magnitude > maxmagnitude) {
					maxmagnitude = magnitude;
					maxMagnatudeEvenIdx = peakidx;

				}
				PeakTmpInfo pti = new PeakTmpInfo();
				pti.setIdx(peakidx);
				pti.setMagnatude((float) magnitude);
				pti.setPeakpos((float) p.getU());
				evenpeaks.add(pti);

				System.out.println("evenpeak\t" + p.getU() + "\t" + magnitude
						+ "\t" + peakidx);

			} else {
				System.out.println("odd\t" + p.getU() + "\t" + p.getR() + "\t"
						+ peakidx);
			}

		}

		// check position of maxEven Idx
		// change if two more lower even peak exist and that have
		// close magnatude
		// code for max even peak is hexaploid case
		// very rare case
		// trim
		double maxmag = 0;
		PeakTmpInfo ptimax = null;
		int cnt = 0;
		for (PeakTmpInfo pti : evenpeaks) {
			//
			if (pti.getIdx() < maxMagnatudeEvenIdx) {
				if (maxmagnitude < (pti.getMagnatude() * 2)) {

					if (pti.getMagnatude() > maxmag) {
						ptimax = pti;
						maxmag = pti.getMagnatude();
					}
					cnt++;
				}
			}
		}
		if (cnt >= 2) {
			//
			maxMagnatudeEvenIdx = ptimax.getIdx();
		}
		// /

		// check largist even peak is based on diploid or tetraploid
		if (maxMagnatudeEvenIdx >= 0) {

			//
			TheoreticalNodematch tnm = new TheoreticalNodematch();
			mmb1 = tnm.matchPeak(pi, maxMagnatudeEvenIdx, map, cnvcount);
			//
			int ploidy = mmb1.getPloidyflg();
			int tp = mmb1.getBestmme().getPurity();
			double tpd = tp * 0.01;
			//
			dataset.setTumorratioFiitiingMatrix(tpd);

		} else {
			// no even number peak was observed, not realistic

		}
		if (pi.getPeaklist().size() <= 1) {
			// no peak detected, analysis was not valid
			dataset.setTumorratioFiitiingMatrix(0f);
		}

		// set copy number
		float ploidy = 0;
		double totallen = 0;
		double totalsum = 0;
		float maxcn = 0;

		//
		for (Peak p : pi.getPeaklist()) {
			for (ChildPeak cp : p.getChildPeaks()) {
				System.out.println(cp.getU() + "\t" + cp.getR() + "\t"
						+ cp.getPeakdist() + "\t" + cp.getTaf() + "\t"
						+ cp.getAaf() + "\t" + cp.getBaf());
			}
		}

		for (List<WaveletIF> wil : dataset.getCapInterval()) {

			for (WaveletIF w : wil) {
				float cn = setCN((CapInterval) w, pi);
				int len = ((CapInterval) w).getLength();
				totallen = totallen + len;
				totalsum = totalsum + (cn * len);
			}
		}
		ploidy = (float) (totalsum / totallen);
		pi.setPloidy(ploidy);

	}

	private float setCN(CapInterval ci, PeaksInfo pi) {

		ci.setCnvtotal(2.0f);
		try {
			List<Peak> pl = pi.getPeaklist();
			Peak peak = pl.get(ci.getPeakIdx());
			return peak.setCN(ci);
		} catch (Exception ex) {

		}
		return 2.0f;
	}

	private Map<SNVHolder, SNVHolderPlusACnv> toMap(AllelicCNV alCNV) {
		List<List<SNVHolderPlusACnv>> list = alCNV.getList();

		Map<SNVHolder, SNVHolderPlusACnv> hsmap = new HashMap<SNVHolder, SNVHolderPlusACnv>();

		for (List<SNVHolderPlusACnv> chlist : list) {
			//
			for (SNVHolderPlusACnv snv : chlist) {
				//
				hsmap.put(snv.getSnv(), snv);

			}

		}
		return hsmap;
	}

}
