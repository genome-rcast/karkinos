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
import java.util.List;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.wavelet.ChildPeak;
import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;

public class TheoreticalNodematch {

	public static void main(String[] arg) {

		//
		TheoreticalNodematch inst = new TheoreticalNodematch();
		List<TheoreticalNodes> list = inst.getPossibleState(
				states_diploid_base, 1.22);
		for (TheoreticalNodes tn : list) {

			System.out.println("ploidy=(" + tn.aAllelePloidy + ":"
					+ tn.bAllelePloidy + ")");
			for (int n = 0; n < 99; n++) {

				float x = tn.imbalance2N[n];
				float y = tn.getDistTo2N()[n];
				System.out.println(x + "\t" + y);

			}

		}

	}

	public static final int diploid = 2;
	public static final int tetraploid = 4;

	public MatchMatrixBean matchPeak(PeaksInfo pi, int maxMagnatudeEvenIdx,
			Map<Integer, PeakAnalysisComponent> map, int cnvcount) {

		//
		Peak pEven = pi.getPeaklist().get(maxMagnatudeEvenIdx);
		PeakAnalysisComponent pacEven = map.get(maxMagnatudeEvenIdx);

		double sdratio = pacEven.sdnormal;
		double sdposition = pEven.getSD();
		List<TheoreticalNodes> diplist = getPossibleState(states_diploid_base,
				pEven.getU());
		List<TheoreticalNodes> tetraplist = getPossibleState(
				states_tetraploid_base, pEven.getU());

		MBResult ret = MatchMatrix.matchMatrix(diplist, tetraplist, pi,
				maxMagnatudeEvenIdx, sdratio, sdposition, map, cnvcount);
		MatchMatrixBean mmb1 = null;
		if (ret.isDubtouse()) {
			mmb1 = ret.getDup();
		} else {
			mmb1 = ret.getTetra();
		}
		mmb1.setDiplist(diplist);
		mmb1.setTetraplist(tetraplist);
		annotatePeak(mmb1, pi, pEven.getU(), maxMagnatudeEvenIdx, sdratio,
				sdposition, map);
		boolean hdtoomach = checkHD(pi);
		if (hdtoomach) {

			//
			if (ret.isDubtouse() && (size(ret.getTetra()) > 0)) {

				mmb1 = ret.getTetra();
				mmb1.setDiplist(diplist);
				mmb1.setTetraplist(tetraplist);
				annotatePeak(mmb1, pi, pEven.getU(), maxMagnatudeEvenIdx,
						sdratio, sdposition, map);

			} else {

				mmb1.takeLargerTP();
				annotatePeak(mmb1, pi, pEven.getU(), maxMagnatudeEvenIdx,
						sdratio, sdposition, map);
			}
		}

		//

		mmb1.setpEvenMax(pEven);
		return mmb1;
		//

		// tetraploid match if even lowpeak or two lower even peak exsit
		// boolean checkTetraPloid = false;
		// int baseploidy = diploid;
		// if(checkTetraPloid){
		//
		// boolean diploid = chackMatch( pi,
		// maxMagnatudeEvenIdx,diplist,tetraplist,sdratio,sdposition);
		// if(!diploid){
		// baseploidy = tetraploid;
		// }
		//
		// }
		// annotate(pEven,baseploidy, pi, maxMagnatudeEvenIdx,
		// map,sdratio,sdposition);
		//

	}

	private int size(MatchMatrixBean tetra) {
		try {
			return tetra.getBestmme().getNodecounter().size();
		} catch (Exception ex) {
		}
		return 0;
	}

	private boolean checkHD(PeaksInfo pi) {

		for (Peak peak : pi.getPeaklist()) {

			//
			for (ChildPeak cp : peak.getChildPeaks()) {
				if (cp.getTaf() == 0) {
					if (cp.getR() > 0.01) {
						return true;
					}
				}
			}

		}
		return false;

	}

	private void annotatePeak(MatchMatrixBean mmb1, PeaksInfo pi, double u,
			int maxMagnatudeEvenIdx, double sdratio, double sdposition,
			Map<Integer, PeakAnalysisComponent> map) {

		int ploidyflg = mmb1.ploidyflg;
		int purity = 1;
		if (mmb1.bestmme != null) {
			purity = mmb1.bestmme.getPurity();
		}
		int[][] states = states_dip;
		if (ploidyflg == tetraploid) {
			states = states_teora;
		}

		List<TheoreticalNodes> theolist = getPossibleState(states, u);
		MatchMatrix.anotate(ploidyflg, purity, theolist, pi,
				maxMagnatudeEvenIdx, sdratio, sdposition, map);

	}

	static final int[][] states_dip = new int[][] { { 5, 1 }, { 6, 0 },
			{ 3, 2 }, { 4, 1 }, { 5, 0 }, { 2, 2 }, { 3, 1 }, { 4, 0 },
			{ 2, 1 }, { 3, 0 }, { 1, 1 }, { 2, 0 }, { 1, 0 }, { 0, 0 } };

	static final int[][] states_teora = new int[][] { { 5, 5 }, { 6, 4 },
			{ 7, 3 }, { 8, 2 }, { 9, 1 }, { 10, 0 }, { 5, 4 }, { 6, 3 },
			{ 7, 2 }, { 8, 1 }, { 9, 0 }, { 4, 4 }, { 5, 3 }, { 6, 2 },
			{ 7, 1 }, { 8, 0 }, { 4, 3 }, { 5, 2 }, { 6, 1 }, { 7, 0 },
			{ 3, 3 }, { 4, 2 }, { 5, 1 }, { 6, 0 }, { 3, 2 }, { 4, 1 },
			{ 5, 0 }, { 2, 2 }, { 3, 1 }, { 4, 0 }, { 2, 1 }, { 3, 0 },
			{ 1, 1 }, { 2, 0 }, { 1, 0 }, { 0, 0 } };

	static final int[][] states_diploid_base = new int[][] { { 3, 1 },
			{ 2, 1 }, { 1, 1 }, { 2, 0 }, { 1, 0 } };

	static final int[][] states_tetraploid_base = new int[][] { { 4, 2 },
			{ 3, 2 }, { 2, 2 }, { 2, 1 }, { 1, 1 }, { 2, 0 } };

	private List<TheoreticalNodes> getPossibleState(int[][] states,
			double peakpos) {

		List<TheoreticalNodes> list = new ArrayList<TheoreticalNodes>();
		for (int[] ploidy : states) {

			TheoreticalNodes tn = new TheoreticalNodes(ploidy, peakpos);
			list.add(tn);

		}
		return list;
	}

}
