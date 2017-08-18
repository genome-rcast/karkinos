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

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import jp.ac.utokyo.rcast.karkinos.wavelet.ChildPeak;
import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;

public class MatchMatrix {

	final static int diploid = 2;
	final static int tetraploid = 4;
	final static int hexaploid = 6;

//	public static MatchMatrixBean matchMatrixHex(List<TheoreticalNodes> dilist,
//			List<TheoreticalNodes> tetraplist, PeaksInfo pi,
//			int maxMagnatudeEvenIdx, double sdratio, double sdposition,
//			Map<Integer, PeakAnalysisComponent> map) {
//
//		// extract vector
//		List<PeakPoint> list = getV(pi, maxMagnatudeEvenIdx, map);
//
//		// ristrict by peak dist
//		List<Float> peakDists = getPeakDist(pi, maxMagnatudeEvenIdx,
//				sdposition, map);
//
//		// extratc theoretical vectors for each diploid, triploid
//		MatchMatrixBean mmbTetra = match(dilist, list, peakDists, sdratio,
//				sdposition, tetraploid,false);
//		MatchMatrixBean mmbHex = match(tetraplist, list, peakDists, sdratio,
//				sdposition, hexaploid,false);
//
//		boolean hextouse = compHex(mmbTetra, mmbHex);
//
//		if (hextouse) {
//
//			mmbTetra.setList(list);
//			return mmbTetra;
//		} else {
//			mmbHex.setList(list);
//			return mmbHex;
//		}
//
//	}

	private static boolean compHex(MatchMatrixBean mmbTetra,
			MatchMatrixBean mmbHex) {

		MatchMatrixEach tetra = mmbTetra.getBestmme();
		MatchMatrixEach hex = mmbHex.getBestmme();

		// for triploidhit requre lower 2n or odd 2 hit
		if (hex.getSumratio() > tetra.getSumratio()) {

			//
			boolean two2 = hex.nodecounter.containsKey("2.0-2.0");
			boolean one1 = hex.nodecounter.containsKey("1.0-1.0");
			if (one1 && two2) {
				return false;
			}

		}

		return true;

	}

	public static MBResult matchMatrix(List<TheoreticalNodes> dilist,
			List<TheoreticalNodes> tetraplist, PeaksInfo pi,
			int maxMagnatudeEvenIdx, double sdratio, double sdposition,
			Map<Integer, PeakAnalysisComponent> map, int cnvcount) {

		MBResult ret = new MBResult();
		
		// extract vector
		List<PeakPoint> list = getV(pi, maxMagnatudeEvenIdx, map);

		// ristrict by peak dist
		List<Float> peakDists = getPeakDist(pi, maxMagnatudeEvenIdx,
				sdposition, map);

		boolean noisy = (cnvcount>250)||(pi.getPeaklist().size()>15);
		boolean notManyPeak = (pi.getPeaklist().size() <=6);
		
		// extratc theoretical vectors for each diploid, triploid
		MatchMatrixBean mmbDi = match(dilist, list, peakDists, sdratio,
				sdposition, diploid,noisy);
		MatchMatrixBean mmbTrip = match(tetraplist, list, peakDists, sdratio,
				sdposition, tetraploid,noisy);
		boolean lowtp = mmbTrip.getBestmme().getPurity()<=18;
		
		boolean diploidtouse = true;
		if(noisy || lowtp || notManyPeak){
			diploidtouse = true;
		}else{
			diploidtouse = comp(mmbDi, mmbTrip);
		}		
		
		mmbDi.setList(list);
		mmbTrip.setList(list);

		ret.setDubtouse(diploidtouse);
		ret.setDup(mmbDi);
		ret.setTetra(mmbTrip);
		
		return ret;

		// PeakAnalysisComponent pac = map.get(maxMagnatudeEvenIdx);
		// if(pac.somaticSNV.getSs().getN()> 50){
		// double meanAFsomatic = pac.somaticSNV.getSs().getMean();
		// double sd = pac.somaticSNV.getSs().getStandardDeviation();
		// double min = meanAFsomatic-sd;
		// double max = meanAFsomatic+sd;
		// if(min<0)min=0;
		// if(max>1)max=1;
		// // restrict by tumor purity from somatic mutation
		// // if enough # of SNV exsist
		// float mintp = (float) (min*2);
		// float maxtp = (float) (max*2);
		// System.out.println(mintp+"\t"+maxtp);
		//
		// }

	}

	private static boolean comp(MatchMatrixBean mmbDi, MatchMatrixBean mmbTrip) {

		MatchMatrixEach di = mmbDi.getBestmme();
		MatchMatrixEach trip = mmbTrip.getBestmme();

		if (di == null && trip == null)
			return true;
		if (di != null && trip == null)
			return true;
		if (di == null && trip != null)
			return false;

		if (trip.getNumhit() <= 2)
			return true;

		// for triploidhit requre lower 2n or odd 2 hit
		if ((trip.getNumhit() >= di.getNumhit())
				&& (trip.getSumratio() > di.getSumratio())) {
						
			
			//
			boolean one1 = trip.nodecounter.containsKey("1.0-1.0");
			boolean two0 = trip.nodecounter.containsKey("2.0-0.0");
			boolean two1 = trip.nodecounter.containsKey("2.0-1.0");
			boolean tetracond = (one1&&two0) || (one1&&two1) || (two0 && two1);
			return !tetracond;

		}

		return true;

	}

	private static MatchMatrixBean match(List<TheoreticalNodes> dilist,
			List<PeakPoint> list, List<Float> peakDists, double sdratio,
			double sdposition, int diploid, boolean noisy) {

		MatchMatrixBean mmb = new MatchMatrixBean();
		mmb.setPloidyflg(diploid);

		List<MatchMatrixEach> macthlist = new ArrayList<MatchMatrixEach>();

		int maxnodehit = 0;
		// purity
		for (int n = 5; n < 100; n++) {

			Map<String, double[]> nodecounter = new LinkedHashMap<String, double[]>();

			double sumdist = 0;
			double sumratio = 0;
			boolean hit = false;
			int numhit = 0;
			int idx = -1;
			boolean largistpeakhit = false;
			for (PeakPoint pp : list) {

				idx++;
				// do not match less thna 1%
				if (pp.peakmagnitude < 0.01) {
					continue;
				}
				//

				TheoreticalNodes tn = closestNode(pp, dilist, n, diploid);
				if (tn == null)
					continue;
				float dist = 0;
				float imblance = 0;
				if (diploid == 2) {
					//
					dist = tn.distTo2N[n];
					imblance = tn.imbalance2N[n];

					// homozigous delation does not exceed 3%
					if ((tn.aAllelePloidy == 0) && (tn.bAllelePloidy == 0)) {
						if (pp.peakmagnitude > 0.03) {
							continue;
						}
					}

				} else {
					//
					dist = tn.distTo4N[n];
					imblance = tn.imbalance4N[n];
				}
				if (outside(pp.peakpos, dist, peakDists, sdposition * 3)) {
					continue;
				}
				if (outside(pp.imbalanceratio, imblance, sdratio * 3)) {
					continue;
				}

				String nodeid = tn.getID();
				if (!nodecounter.containsKey(nodeid)) {

					// hit
					hit = true;
					if(idx==0){
						largistpeakhit = true;
					}
					double eucliddist = dist(pp.peakpos - dist,
							pp.imbalanceratio - imblance);
					sumdist = sumdist + eucliddist;
					sumratio = sumratio + pp.peakmagnitude;
					numhit++;

					nodecounter.put(nodeid, new double[] { eucliddist,
							pp.peakmagnitude });

				}

			}
			if (hit) {
				MatchMatrixEach mme = new MatchMatrixEach();
				mme.setNodecounter(nodecounter);
				mme.setPurity(n);
				mme.setSumdist(sumdist);
				mme.setSumratio(sumratio);
				macthlist.add(mme);
				mme.setNumhit(numhit);
				if (nodecounter.size() > maxnodehit) {
					maxnodehit = nodecounter.size();
				}
			}

		}

		boolean noisydiploid = (diploid==2)&&noisy;
		//
		if (maxnodehit == 0) {

			MatchMatrixEach mme = new MatchMatrixEach();
			mme.setPurity(20);
			macthlist.add(mme);
			mme.setNumhit(0);
			mmb.setBestmme(mme);
			return mmb;

		} else if (maxnodehit == 1 || noisydiploid) {

			// if hit = 1
			// take loh cand if exsist
			//
			MatchMatrixEach bestmme = getBestLOH(macthlist);
			mmb.setMacthlist(macthlist);
			mmb.setBestmme(bestmme);
			//
			return mmb;
		}

		Map<Long, List<MatchMatrixEach>> map = new LinkedHashMap<Long, List<MatchMatrixEach>>();

		for (MatchMatrixEach mme : macthlist) {

			//
			if (mme.getNodecounter().size() == maxnodehit) {

				long sr = Math.round(mme.sumratio * 100);
				List<MatchMatrixEach> llist = null;
				if (map.containsKey(sr)) {
					llist = map.get(sr);
				} else {
					llist = new ArrayList<MatchMatrixEach>();
					map.put(sr, llist);
				}
				//
				llist.add(mme);

			}

		}

		List<MatchMatrixEach> bestmmeList = new ArrayList<MatchMatrixEach>();
		// extract max sumratio candidate among max hit vector
		Iterator<Long> ite = map.keySet().iterator();
		while (ite.hasNext()) {
			//
			Long key = ite.next();
			List<MatchMatrixEach> llist2 = map.get(key);
			MatchMatrixEach bestmme = null;
			double mindist = Double.MAX_VALUE;
			for (MatchMatrixEach mme : llist2) {

				//
				if (mme.sumdist < mindist) {
					mindist = mme.sumdist;
					bestmme = mme;
				}

			}
			bestmmeList.add(bestmme);

		}

		// chose best hit among different tumor purity cand
		MatchMatrixEach maxpurityMme = null;
		MatchMatrixEach maxratioMme = null;
		// take max tumor purity cand
		int maxpurity = 0;
		double maxratio = 0;
		for (MatchMatrixEach mme : bestmmeList) {

			//
			if (mme.purity > maxpurity) {

				//
				maxpurity = mme.purity;
				maxpurityMme = mme;
			}

			if (mme.sumratio > maxratio) {

				//
				maxratio = mme.sumratio;
				maxratioMme = mme;
			}

		}

		MatchMatrixEach bestmme = null;
		if (maxpurityMme == maxratioMme) {
			bestmme = maxpurityMme;
		} else {

			bestmme = compare(maxpurityMme, maxratioMme);

		}
		//
		mmb.setMacthlist(macthlist);
		mmb.setBestmme(bestmme);
		//
		return mmb;
	}

	private static MatchMatrixEach getBestLOH(List<MatchMatrixEach> macthlist) {

		boolean containloh = false;
		List<MatchMatrixEach> lohlist = new ArrayList<MatchMatrixEach>();
		for (MatchMatrixEach mme : macthlist) {
			if (mme.nodecounter.containsKey("1.0-0.0")) {
				lohlist.add(mme);
				containloh = true;
			}
		}
		if (containloh) {
			return bestnode2(lohlist);
		} else {
			return bestnode2(macthlist);
		}

	}

	private static MatchMatrixEach bestnode(List<MatchMatrixEach> list) {

		double dist = Double.MAX_VALUE;
		MatchMatrixEach ret = null;
		for (MatchMatrixEach mme : list) {
			if (mme.sumdist < dist) {
				dist = mme.sumdist;
				ret = mme;
			}
		}
		return ret;
	}

	private static MatchMatrixEach bestnode2(List<MatchMatrixEach> list) {

		Map<String, List<MatchMatrixEach>> map = new LinkedHashMap<String, List<MatchMatrixEach>>();
		int maxpurity = 0;
		String maxpurityKey = null;
		for (MatchMatrixEach mme : list) {

			Iterator<String> ite = mme.nodecounter.keySet().iterator();
			if (ite.hasNext()) {

				String key = ite.next();

				if (mme.getPurity() > maxpurity) {
					maxpurity = mme.getPurity();
					maxpurityKey = key;
				}

				List<MatchMatrixEach> listr = null;
				if (map.containsKey(key)) {
					listr = map.get(key);
				} else {
					listr = new ArrayList<MatchMatrixEach>();
					map.put(key, listr);
				}
				listr.add(mme);
			}

		}
		//
		return bestnode(map.get(maxpurityKey));

	}

	private static MatchMatrixEach compare(MatchMatrixEach maxpurityMme,
			MatchMatrixEach maxratioMme) {

		double ratio0 = maxpurityMme.getSumratio();
		double ratio1 = maxratioMme.getSumratio();
		if (ratio0 == 0)
			return maxpurityMme;// shondnt be, just for safety
		// if ratio is too small for max purity candidate
		// return max ratio cand
		boolean ratio0toosamall = (ratio0 / ratio1) < 0.01;
		if (ratio0toosamall) {
			if (maxratioMme != null) {
				return maxratioMme;
			}
		}
		return maxpurityMme;
	}

	private static TheoreticalNodes closestNode(PeakPoint pp,
			List<TheoreticalNodes> dilist, int n, int diploid) {

		double mindist = Double.MAX_VALUE;
		TheoreticalNodes tnbest = null;
		for (TheoreticalNodes tn : dilist) {
			//
			float dist = 0;
			float imblance = 0;
			if (diploid == 2) {
				//
				dist = tn.distTo2N[n];
				imblance = tn.imbalance2N[n];
			} else {
				//
				dist = tn.distTo4N[n];
				imblance = tn.imbalance4N[n];
			}
			double eucliddist = dist(pp.peakpos - dist, pp.imbalanceratio
					- imblance);
			if (mindist > eucliddist) {
				tnbest = tn;
				mindist = eucliddist;
			}
		}
		return tnbest;
	}

	private static double dist(float f, float g) {

		return Math.sqrt(pow2(f) + pow2(g));
	}

	private static double pow2(float g) {

		return Math.pow(g, 2);
	}

	private static boolean outside(float imblance, float imbalanceratio,
			double sdratio) {

		return Math.abs(imblance - imbalanceratio) > (3 * sdratio);

	}

	private static boolean outside(float dist, float distth,
			List<Float> peakDists, double sdposition) {
		int[] nn = new int[] { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7, 7,
				-8, 8, -9, 9, -10, 10 };
		for (int n : nn) {
			for (Float f : peakDists) {
				//
				if (Math.abs(distth - (f * n)) < sdposition) {

					if (Math.abs(dist - distth) < sdposition) {
						return false;
					}
				}

			}
		}
		return true;
	}

	private static List<PeakPoint> getV(PeaksInfo pi, int maxMagnatudeEvenIdx,
			Map<Integer, PeakAnalysisComponent> map) {

		List<PeakPoint> list = new ArrayList<PeakPoint>();
		Peak pEven = pi.getPeaklist().get(maxMagnatudeEvenIdx);
		PeakAnalysisComponent pacEven = map.get(maxMagnatudeEvenIdx);

		for (int n = 0; n < pi.getPeaklist().size(); n++) {

//			if (n == maxMagnatudeEvenIdx) {
//				continue;
//			}
			boolean peakself = (n == maxMagnatudeEvenIdx);
			Peak peak = pi.getPeaklist().get(n);
			if (peak.isArtifitial())
				continue;
			PeakAnalysisComponent pac = map.get(n);
			if (pac == null)
				continue;

			//
			if (pac.complexpeak) {


				PeakPoint pp2 = new PeakPoint();
				float peakpos2 = (float) (peak.getU() - pEven.getU());
				float imbalanceratio2 = (float) (pac.oddpd - pacEven.hpeakdistance);
				float peakmagnitude2 = (float) (peak.getR() * pac.oddr);
				pp2.setPeakpos(peakpos2);
				pp2.setImbalanceratio(imbalanceratio2);
				pp2.setPeakmagnitude(peakmagnitude2);
				
				if (peakself) {
					
					if(imbalanceratio2>0.25){
						list.add(pp2);
					}
					continue;
					
				}else{
					list.add(pp2);
				}
							
				PeakPoint pp1 = new PeakPoint();
				float peakpos1 = (float) (peak.getU() - pEven.getU());
				float imbalanceratio1 = (float) (pac.evenpd - pacEven.hpeakdistance);
				float peakmagnitude1 = (float) (peak.getR() * pac.evenr);
				pp1.setPeakpos(peakpos1);
				pp1.setImbalanceratio(imbalanceratio1);
				pp1.setPeakmagnitude(peakmagnitude1);
				list.add(pp1);
				

			} else {
				
				if (peakself) {
					continue;
				}
				
				PeakPoint pp = new PeakPoint();
				float peakpos = (float) (peak.getU() - pEven.getU());
				float imbalanceratio = (float) (pac.hpeakdistance - pacEven.hpeakdistance);
				float peakmagnitude = (float) (peak.getR());
				pp.setPeakpos(peakpos);
				pp.setImbalanceratio(imbalanceratio);
				pp.setPeakmagnitude(peakmagnitude);
				list.add(pp);

			}

		}
		sort(list);
		list = top10(list);
		return list;
	}

	private static List<PeakPoint> top10(List<PeakPoint> list) {
		
		List<PeakPoint> rlist = new ArrayList<PeakPoint>();
		if(list.size()<=10){
			return list;
		}else{
			int i = 0;
			while(i<10){
				rlist.add(list.get(i));
				i++;
			}
		}
		return rlist;
	}

	private static void sort(List<PeakPoint> list) {

		Collections.sort(list);

	}

	private static List<Float> getPeakDist(PeaksInfo pi,
			int maxMagnatudeEvenIdx, double sd,
			Map<Integer, PeakAnalysisComponent> map) {

		//
		List<Float> pdist = new ArrayList<Float>();
		Peak pEven = pi.getPeaklist().get(maxMagnatudeEvenIdx);
		int[] nn = new int[] { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6, -7, 7,
				-8, 8, -9, 9, -10, 10 };
		for (int n : nn) {

			int idx = maxMagnatudeEvenIdx + n;
			if (idx < 0)
				continue;
			if (idx >= pi.getPeaklist().size())
				continue;

			//
			Peak peak = pi.getPeaklist().get(idx);
			float dist = (float) Math.abs(pEven.getU() - peak.getU());
			if (pdist.size() == 0) {
				pdist.add(dist);
			} else {

				float dist2 = notinmultiple((float) dist, pdist, sd);
				if (dist2 < 0) {
					pdist.add(dist2);
				}

			}

		}
		return pdist;

	}

	private static float notinmultiple(float f, List<Float> pdist2, double sd) {

		for (Float fd : pdist2) {

			for (int n = 1; n <= 10; n++) {

				//
				if (include(fd, f / n, sd)) {
					return Math.abs((f / n));
				}

			}

		}
		return -100;

	}

	private static boolean include(float f, float f2, double sd) {

		return Math.abs(f - f2) < sd;

	}

	public static void anotate(int baseploidy, int purity,
			List<TheoreticalNodes> list, PeaksInfo pi, int maxMagnatudeEvenIdx,
			double sdratio, double sdposition,
			Map<Integer, PeakAnalysisComponent> map) {

		//
		int idx = 0;
		for (Peak peak : pi.getPeaklist()) {
			//
			peak.setPac(map.get(idx));
			idx++;
		}

		//
		idx = 0;
		Peak maxEven = pi.getPeaklist().get(maxMagnatudeEvenIdx);
		for (Peak peak : pi.getPeaklist()) {

			if (peak == maxEven) {
				matchOriginNode(peak, purity, baseploidy);
				continue;
			}
			if (peak.getPac() == null) {
				continue;
			}

			for (ChildPeak cp : peak.getChildPeaks()) {

				match(maxEven, cp, list, baseploidy, purity,
						maxMagnatudeEvenIdx, idx, sdratio, sdposition);

			}
			idx++;
		}

	}

	private static void matchOriginNode(Peak peak, int purity, int baseploidy) {

		double purityd = purity * 0.01;
		for (ChildPeak cp : peak.getChildPeaks()) {

			if (baseploidy == diploid) {
				cp.setAaf(1);
				cp.setBaf(1);
				cp.setTaf(2);

				//
				double r = (1 + purityd) / 2;
				double margin = 0.1;
				double d1 = Math.abs((cp.getPeakdist() - r));
				double d2 = Math.abs((cp.getPeakdist() - (0.5 + margin)));

				//
				if ((cp.getPeakdist() > 0.65) && (d1 < d2)) {
					if (cp.isNotSexChrom()) {

						cp.setAaf(2);
						cp.setBaf(0);
						cp.setTaf(2);
					}

				}

			} else {

				cp.setAaf(2);
				cp.setBaf(2);
				cp.setTaf(4);

				double r = (2 * purityd + 1) / (2 * purityd + 2);
				double margin = 0.1;
				double d1 = Math.abs((cp.getPeakdist() - r));
				double d2 = Math.abs((cp.getPeakdist() - (0.5 + margin)));

				//
				if ((cp.getPeakdist() > 0.65) && (d1 < d2)) {
					if (cp.isNotSexChrom()) {
						cp.setAaf(3);
						cp.setBaf(1);
						cp.setTaf(4);
					}

				}
			}
		}

	}

	private static void match(Peak maxEven, ChildPeak cp,
			List<TheoreticalNodes> list, int baseploidy, int purity,
			int maxMagnatudeEvenIdx, int idx, double sdratio, double sdposition) {

		TheoreticalNodes closest = getClosest(maxEven, cp, baseploidy, purity,
				list, sdratio, sdposition);

		//
		// hit to thoretical
		cp.setAaf(closest.getaAllelePloidy());
		cp.setBaf(closest.getbAllelePloidy());
		cp.setTaf(closest.getTotalploidy());

	}

	private static TheoreticalNodes getClosest(Peak maxEven, ChildPeak cp,
			int baseploidy, int purity, List<TheoreticalNodes> list,
			double sdratio, double sdposition) {

		float af0 = (float) cp.getPeakdist() - maxEven.getPac().hpeakdistance;
		float u0 = (float) (cp.getU() - maxEven.getU());

		int idx = 0;
		int minidx = 0;
		float mindist = Float.MAX_VALUE;
		boolean hit = false;

		double unitpeakdist = 0.5;
		for (TheoreticalNodes tn : list) {

			float af = 0f;
			float u = 0f;
			//
			if (baseploidy == diploid) {
				af = tn.imbalance2N[purity];
				u = tn.distTo2N[purity];
				if (tn.getTotalploidy() == 1) {
					unitpeakdist = Math.abs(u);
				}

			} else {
				af = tn.imbalance4N[purity];
				u = tn.distTo4N[purity];
				if (tn.getTotalploidy() == 3) {
					unitpeakdist = Math.abs(u);
				}
			}
			boolean match1 = Math.abs(u0 - u) < (sdposition * 3);
			if (!match1) {
				idx++;
				continue;
			}
			//
			double eucliddist = dist(u0 - u, af0 - af);
			if (eucliddist < mindist) {
				mindist = (float) eucliddist;
				minidx = idx;
				tn.setDiffu(u0 - u);
				tn.setDiffpd(af0 - af);
				hit = true;
			}
			//
			idx++;

		}
		if (hit) {
			TheoreticalNodes tnr = list.get(minidx);
			return tnr;
		} else {

			TheoreticalNodes notartifitial = new TheoreticalNodes();
			double totalploidy = 0;
			if (baseploidy == diploid) {
				totalploidy = 2 + (u0 / unitpeakdist);
			} else {
				totalploidy = 4 + (u0 / unitpeakdist);
			}
			notartifitial.setTotalploidy((float) totalploidy);
			double afAdj = (double) ((double) (cp.getPeakdist() - 0.5) / (double) (purity * 0.01));
			if (afAdj > 0.5)
				afAdj = 0.5;
			double aaf = (afAdj + 0.5) * totalploidy;
			double baf = totalploidy - aaf;
			notartifitial.setaAllelePloidy((float) aaf);
			notartifitial.setbAllelePloidy((float) baf);
			return notartifitial;

		}

	}
}
