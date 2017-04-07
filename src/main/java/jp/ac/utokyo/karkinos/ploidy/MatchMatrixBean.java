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

import java.util.List;

import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;

public class MatchMatrixBean {

	int ploidyflg;
	List<MatchMatrixEach> macthlist;
	MatchMatrixEach bestmme;

	Peak pEvenMax;

	public Peak getpEvenMax() {
		return pEvenMax;
	}

	public void setpEvenMax(Peak pEvenMax) {
		this.pEvenMax = pEvenMax;
	}

	// for display
	List<PeakPoint> list;

	public List<TheoreticalNodes> getDiplist() {
		return diplist;
	}

	public void setDiplist(List<TheoreticalNodes> diplist) {
		this.diplist = diplist;
	}

	public List<TheoreticalNodes> getTetraplist() {
		return tetraplist;
	}

	public void setTetraplist(List<TheoreticalNodes> tetraplist) {
		this.tetraplist = tetraplist;
	}

	List<TheoreticalNodes> diplist = null;
	List<TheoreticalNodes> tetraplist = null;

	public List<PeakPoint> getList() {
		return list;
	}

	public void setList(List<PeakPoint> list) {
		this.list = list;
	}

	public MatchMatrixEach getBestmme() {
		return bestmme;
	}

	public void setBestmme(MatchMatrixEach bestmme) {
		this.bestmme = bestmme;
	}

	public List<MatchMatrixEach> getMacthlist() {
		return macthlist;
	}

	public void setMacthlist(List<MatchMatrixEach> macthlist) {
		this.macthlist = macthlist;
	}

	public int getPloidyflg() {
		return ploidyflg;
	}

	public void setPloidyflg(int ploidyflg) {
		this.ploidyflg = ploidyflg;
	}

	int nummatch = 0;
	double summatchdist = 0;

	public int getNummatch() {
		return nummatch;
	}

	public void setNummatch(int nummatch) {
		this.nummatch = nummatch;
	}

	public double getSummatchdist() {
		return summatchdist;
	}

	public void setSummatchdist(double summatchdist) {
		this.summatchdist = summatchdist;
	}

	public double getInterval() {
		if (pEvenMax != null) {

			double pp = pEvenMax.getU();
			double pt = 0.01 * getBestmme().getPurity();
			return getUnitDist(pt, getPloidyflg(), pp);
		}
		return 0;
	}

	private static float getUnitDist(double pt, int baseploidy, double peakpos) {

		if (baseploidy == 2) {
			return (float) ((pt / 2) * peakpos);
		} else {
			// y = -0.15x2+0.4x
			return (float) ((-0.15 * Math.pow(pt, 2) + (0.4 * pt)) * peakpos);
		}

	}

	public void takeLargerTP() {
		int purity = 0;
		if (macthlist != null) {
			for (MatchMatrixEach mme : macthlist) {
				if(mme.purity>purity){
					purity = mme.purity;
					bestmme = mme;
				}
				
			}
		}

	}
}
