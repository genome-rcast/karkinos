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
import java.util.Collections;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.wavelet.FunctionRegression;

public class AllelicCNV {

	DataSet dataset;
	List<List<SNVHolderPlusACnv>> list = new ArrayList<List<SNVHolderPlusACnv>>();
	ReadsSummary readsSummary;

	public AllelicCNV(DataSet dataset, ReadsSummary readsSummary) {

		this.dataset = dataset;
		this.readsSummary = readsSummary;
		analysis(dataset);
		setIntervalList();

	}

	List<CopyNumberInterval> allelicLOHhigh;
	List<CopyNumberInterval> allelicLOHLow;
	List<CopyNumberInterval> hdList;

	private void setIntervalList() {

		allelicLOHLow = getLOH(false);
		allelicLOHhigh = getLOH(true);
		hdList = analyse(allelicLOHLow, allelicLOHhigh);
		allelicLOHLow = takeLowOnly(allelicLOHLow);
		allelicLOHhigh = takeHighOnly(allelicLOHhigh);
	}

	public List<CopyNumberInterval> getSortedList() {
		List<CopyNumberInterval> list = new ArrayList<CopyNumberInterval>();

		list.addAll(allelicLOHLow);
		list.addAll(allelicLOHhigh);
		list.addAll(hdList);
		sort(list);
		return list;
	}

	private void sort(List<CopyNumberInterval> list) {

		Collections.sort(list, new MyComp());

	}

	private List<CopyNumberInterval> takeHighOnly(
			List<CopyNumberInterval> allelicCNLIST) {

		List<CopyNumberInterval> list = new ArrayList<CopyNumberInterval>();
		for (CopyNumberInterval cni : allelicCNLIST) {
			if (cni.getCopynumber() > 2) {
				list.add(cni);
			}
		}
		return list;
	}

	private List<CopyNumberInterval> takeLowOnly(
			List<CopyNumberInterval> allelicCNLIST) {
		List<CopyNumberInterval> list = new ArrayList<CopyNumberInterval>();
		for (CopyNumberInterval cni : allelicCNLIST) {
			if (cni.getCopynumber() < 2) {
				list.add(cni);
			}
		}
		return list;
	}

	public List<CopyNumberInterval> getAllelicLOHhigh() {
		return allelicLOHhigh;
	}

	public List<CopyNumberInterval> getAllelicLOHLow() {
		return allelicLOHLow;
	}

	public List<CopyNumberInterval> getHdList() {
		return hdList;
	}

	private List<CopyNumberInterval> analyse(
			List<CopyNumberInterval> allelicLOHLow,
			List<CopyNumberInterval> allelicLOHhigh) {

		List<CopyNumberInterval> hdlist = new ArrayList<CopyNumberInterval>();

		for (CopyNumberInterval cni : allelicLOHLow) {
			//
			for (CopyNumberInterval cnihigh : allelicLOHhigh) {
				if (include(cnihigh, cni)) {
					if (cnihigh.getCopynumber() == 1) {
						cnihigh.setHdelation(true);
						hdlist.add(cnihigh);
						cnihigh.setCopynumber(0);
					}
				}
			}
		}
		return hdlist;
	}

	private boolean include(CopyNumberInterval cnihigh, CopyNumberInterval cni) {

		if (!cnihigh.getChr().equals(cni.getChr())) {
			return false;
		}
		if (cnihigh.getCopynumber() != cni.getCopynumber()) {
			return false;
		}
		if ((cni.getStart() <= cnihigh.getStart())
				&& (cni.getEnd() >= cnihigh.getEnd())) {

			if (len(cni) > len(cnihigh)) {

				return true;
			}
		}
		return false;
	}

	private int len(CopyNumberInterval cni) {

		return Math.abs(cni.getEnd() - cni.getStart());

	}

	public List<CopyNumberInterval> getLOH(boolean higher) {

		List<CopyNumberInterval> copyNumberIntervalList = new ArrayList<CopyNumberInterval>();

		//
		for (List<SNVHolderPlusACnv> plist : list) {
			int statusb4 = -1;
			CopyNumberInterval cni = null;
			int posb4 = 0;
			int noSNP = 0;
			for (SNVHolderPlusACnv sc : plist) {
				//
				int hmmv = 0;
				if (higher) {
					hmmv = (int) sc.highera.hmmval;
				} else {
					hmmv = (int) sc.lowera.hmmval;
				}
				if (hmmv != statusb4) {

					if (cni != null && cni.getCopynumber() != 2) {
						cni.setEnd(posb4);
						cni.setNoSNP(noSNP);
						copyNumberIntervalList.add(cni);
					}
					cni = new CopyNumberInterval();
					cni.setChr(sc.getSnv().getChr());
					cni.setStart(sc.getSnv().getPos());
					cni.setCopynumber(hmmv);
					cni.setAllelic(true);
					noSNP = 0;
				}
				noSNP++;
				statusb4 = hmmv;
				posb4 = sc.getSnv().getPos();
			}
			if (cni != null && cni.getCopynumber() != 2) {

				if (posb4 > cni.getStart()) {
					cni.setEnd(posb4);
					copyNumberIntervalList.add(cni);
				}
			}
		}
		//

		return copyNumberIntervalList;

	}

	private void analysis(DataSet dataset) {

		// set low data and GC adjusted
		// FunctionRegression fr = dataset.getFunctionRegression();
		//
		int minbaitlen = 1000;
		List<double[]> xyzListHigh = new ArrayList<double[]>();
		List<double[]> xyzListLow = new ArrayList<double[]>();

		double ntratio = 1;
		try {

			int nrc = readsSummary.getNormalCounter().getTotalOnTarget();
			int trc = readsSummary.getTumorCounter().getTotalOnTarget();
			ntratio = trc / nrc;
			if (Double.isInfinite(ntratio) || Double.isNaN(ntratio)
					|| (ntratio < 0.2) || (ntratio > 5)) {
				ntratio = 1;
			}

		} catch (Exception ex) {

		}

		for (SNVHolder snv : dataset.getSnvlist()) {

			if (snv.isHetroSNP()) {

				SNVHolderPlusACnv sp = new SNVHolderPlusACnv(snv,ntratio);
				if (sp.valid) {
					CapInterval civ = snv.getCi();
					double x = civ.getLength();
					double y = civ.getCgParcent();
					double[] arLow = new double[3];
					arLow[0] = x;
					arLow[1] = y;
					arLow[2] = sp.lowera.row;
					xyzListLow.add(arLow);

					double[] arHigh = new double[3];
					arHigh[0] = x;
					arHigh[1] = y;
					arHigh[2] = sp.highera.row;
					xyzListHigh.add(arHigh);

					if (civ.getLength() < minbaitlen) {
						minbaitlen = civ.getLength();
					}
				}
			}

		}
		// getGCAdjstFunction
		System.out.println("GC AL High");
		FunctionRegression frHigh = new FunctionRegression(xyzListHigh,
				minbaitlen);
		System.out.println("GC AL Low");
		FunctionRegression frLow = new FunctionRegression(xyzListLow,
				minbaitlen);

		List<SNVHolderPlusACnv> sublist = null;
		String chrom = "";
		for (SNVHolder snv : dataset.getSnvlist()) {

			if (snv.isHetroSNP()) {

				if (!snv.getChr().equals(chrom)) {

					chrom = snv.getChr();
					if (sublist != null) {
						list.add(sublist);
					}
					sublist = new ArrayList<SNVHolderPlusACnv>();
				}

				SNVHolderPlusACnv sp = new SNVHolderPlusACnv(snv,ntratio);
				CapInterval civ = snv.getCi();
				double x = civ.getLength();
				double y = civ.getCgParcent();
				//

				if (sp.valid) {
					sp.lowera.gcadjusted = (float) frLow.getAdjustedZ(x, y,
							sp.lowera.row);
					sp.highera.gcadjusted = (float) frHigh.getAdjustedZ(x, y,
							sp.highera.row);
					sublist.add(sp);
				}
				//

			}
		}
		if (sublist != null) {
			list.add(sublist);
		}

		// //
		// wavelet transform
		try {
			WaveletDenoizeACNV.calc(list);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// HMM anotation
		try {
			HMMACNV.calc(list);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public List<List<SNVHolderPlusACnv>> getList() {
		return list;
	}

}
