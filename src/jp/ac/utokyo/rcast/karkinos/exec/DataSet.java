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
package jp.ac.utokyo.rcast.karkinos.exec;

import htsjdk.samtools.SAMRecord;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.distribution.AnalyseDist;
import jp.ac.utokyo.rcast.karkinos.distribution.DataHolderByCN;
import jp.ac.utokyo.rcast.karkinos.utils.Interval;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;
import jp.ac.utokyo.rcast.karkinos.wavelet.FunctionRegression;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class DataSet implements java.io.Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 6845809656002092214L;

	Map<String, Integer> seqIndex = null;
	boolean useAvearageNormal = false;

	public void setUseAvearageNormal(boolean useAvearageNormal) {
		this.useAvearageNormal = useAvearageNormal;
	}

	public DataSet(Map<String, Integer> _seqIndex) {
		seqIndex = _seqIndex;
	}

	public List<SNVHolder> getSnvlist() {
		return snvlist;
	}

	Coverage normal = new Coverage();
	Coverage tumor = new Coverage();

	public void bailAnalysis() {

		normal.analysisBaitSampling();
		tumor.analysisBaitSampling();
	}

	public Coverage getNormal() {
		return normal;
	}

	public Coverage getTumor() {
		return tumor;
	}

	List<SNVHolder> snvlist = new ArrayList<SNVHolder>();

	public void addSNVHolder(SNVHolder snvholder) {
		snvlist.add(snvholder);
	}

	public void assginCaptureInterval() {

		int assumedReadlen = 100;
		for (SNVHolder holder : snvlist) {

			CapInterval ci = ch
					.getCapInterval(holder.getChr(), holder.getPos());
			if (ci == null) {
				ci = ch.getOverlapping(holder.getChr(), holder.getPos()
						- assumedReadlen, holder.getPos() + assumedReadlen);
			}
			if (ci == null) {
				ci = ch.getOverlapping(holder.getChr(), holder.getPos()
						- KarkinosProp.baitmergin, holder.getPos()
						+ KarkinosProp.baitmergin);
			}
			if (ci == null) {
				ci = ch.getOverlapping(holder.getChr(), holder.getPos()
						- (KarkinosProp.baitmergin * 2), holder.getPos()
						+ (KarkinosProp.baitmergin * 2));
			}

			holder.setCi(ci);

		}
	}

	public void assgindbSNP(DbSNPAnnotation dbSNPAnnotation) {

		for (SNVHolder holder : snvlist) {

			DbSNPBean dbSNPbean = dbSNPAnnotation.lookup(holder.getChr(),
					holder.getPos());
			holder.setDbSNPbean(dbSNPbean);
			boolean hetroSNP = false;
			if (dbSNPbean != null && !dbSNPbean.isCosmic()) {
				hetroSNP = (holder.normal.getRatio() >= KarkinosProp.hetroSNPMin)
						&& (holder.normal.getRatio() <= KarkinosProp.hetroSNPMax);

			} else {

				if (holder.normal.totalcnt >= 20) {
					hetroSNP = (holder.normal.getRatio() >= KarkinosProp.hetroSNPMin)
							&& (holder.normal.getRatio() <= KarkinosProp.hetroSNPMax);
				}

			}
			holder.serHetroSNP(hetroSNP);

		}

	}

	public void setBinReads(Interval iv, int normalcount, int tumorcount) {

		normal.addBinCount(iv, normalcount);
		tumor.addBinCount(iv, tumorcount);

	}

	public OntagetInfo setNomalCoverageInfo(SAMRecord sam) {

		//boolean inTarget = false;
		OntagetInfo ret = null;
		normal.setCoverageInfo(sam);
		//if (targetSeq) {
			ret = normal.setCaptureInfo(sam, ch);
		//}
		return ret;

	}

	public OntagetInfo setTumorCoverageInfo(SAMRecord sam) {

		//boolean inTarget = false;
		OntagetInfo ret = null;
		tumor.setCoverageInfo(sam);
		//if (targetSeq) {
			ret = tumor.setCaptureInfo(sam, ch);
		//}
		return ret;
	}

	public CaptureHolder getCh() {
		return ch;
	}

	// capture region
	CaptureHolder ch = new CaptureHolder();
	boolean targetSeq = false;

	public void loadTargetBed(String targetRegion, TwoBitGenomeReader tgr)
			throws IOException {
		ch.loadTargetBed(targetRegion, tgr);
		if (targetRegion != null) {
			targetSeq = true;
		}
	}

	List<List<WaveletIF>> capinterval = null;
	boolean capintervalinit = false;

	public List<List<WaveletIF>> getCapInterval() {
		if (!capintervalinit) {
			return _getCapInterval();
		} else {
			return capinterval;
		}

	}

	List<CopyNumberInterval> copyNumberIntervalList = null;

	public List<CopyNumberInterval> getCopyNumberIntervalList(int mode) {

		if (copyNumberIntervalList == null) {
			copyNumberIntervalList = calcCopyNumberIntervalList(mode);
		}
		return copyNumberIntervalList;
//		return calcCopyNumberIntervalList(mode);

	}

	public static final int MODE_HMM = 0;
	public static final int MODE_CorrelVaridate = 1;

	private List<CopyNumberInterval> calcCopyNumberIntervalList(int mode) {

		//
		List<CapInterval> alist = new ArrayList<CapInterval>();
		for (List<WaveletIF> list : getCapInterval()) {

			float ab4 = 0;
			float bb4 = 0;
			CapInterval cib4 = null;
			int cnt = 0;
			for (WaveletIF wi : list) {

				boolean isFirst = cnt == 0;
				boolean isLast = cnt == list.size() - 1;
				cnt++;
				CapInterval ci = (CapInterval) wi;
				ci.setStartChrom(isFirst);
				ci.setEndChrom(isLast);

				// double d = ci.getHMMValue();
				// if (mode == MODE_CorrelVaridate) {
				// d = ci.getVaridateVal();
				// }
				// if (mode == 9) {
				// d = ci.getAafreq() * ci.getBafreq();
				// }
				// float n = (float) (d / 0.5);
				// boolean statechage = ((nb4 != n));
				boolean statesame = (ci.getAafreq() == ab4)
						&& (ci.getBafreq() == bb4);
				boolean statechage = !statesame;

			
			
				
				if (statechage || isFirst || isLast) {
					//
					if (cib4 != null && !isLast) {
						alist.add(cib4);
					}
					alist.add(ci);

				}
				ab4 = ci.getAafreq();
				bb4 = ci.getBafreq();
				cib4 = ci;

			}
		}

		List<CopyNumberInterval> copyNumberIntervalList = new ArrayList<CopyNumberInterval>();

		for (int m = 0; m + 1 < alist.size(); m = m + 2) {

			CapInterval cistart = alist.get(m);
			CapInterval ciend = alist.get(m + 1);
			CopyNumberInterval cni = new CopyNumberInterval();
			cni.setChr(cistart.getChr());
			if (cistart.isStartChrom()) {
				cni.setStart(1);
			} else {
				cni.setStart(cistart.getStart());
			}
			if (ciend.isEndChrom()) {
				int end = ciend.getEnd();
				if (seqIndex != null && seqIndex.containsKey(ciend.chr)) {
					cni.setEnd(seqIndex.get(ciend.chr));
				}
				cni.setEnd(end);
			} else {
				cni.setEnd(ciend.getEnd());
			}
			double d = cistart.getHMMValue();
			if (mode == MODE_CorrelVaridate) {
				d = cistart.getVaridateVal();
			}
			if (mode == 9) {
				d = cistart.getAafreq() + cistart.getBafreq();
			}
			//float copynumber = (float) (d / 0.5);
			float copynumber = (float)d;
			cni.setCopynumber(copynumber);
			cni.setAaf(cistart.getAafreq());
			cni.setBaf(cistart.getBafreq());
			boolean contain = false;
			for (CopyNumberInterval cn : copyNumberIntervalList) {
				if (cn.getChr().equals(cni.getChr())) {
					if (cn.getStart() == cni.getStart()) {
						// if(cn.getEnd() == cni.getEnd()){
						if ((cn.aaf == cni.aaf) && (cn.baf == cni.baf)) {
							contain = true;
						}
						// }
					}
				}
			}
			if (contain == false) {
				copyNumberIntervalList.add(cni);
			}
		}

		return copyNumberIntervalList;
	}

	Set<String> chrSet = new LinkedHashSet<String>();

	public List<List<WaveletIF>> _getCapInterval() {

		List<List<WaveletIF>> plist = new ArrayList<List<WaveletIF>>();
		List<WaveletIF> list = new ArrayList<WaveletIF>();
		Map<CapInterval, Coverage.Container> m1 = normal.capregion;
		Map<CapInterval, Coverage.Container> m2 = tumor.capregion;
		Iterator<CapInterval> ite = m1.keySet().iterator();
		System.out.println("normal million adjust=" + normal.millionAdjust());
		System.out.println("tumor million adjust=" + tumor.millionAdjust());

		// /
		SummaryStatistics ss = new SummaryStatistics();
		SummaryStatistics tnss = new SummaryStatistics();
		while (ite.hasNext()) {

			CapInterval key = ite.next();
			Coverage.Container co = m1.get(key);
			Coverage.Container co2 = m2.get(key);
			long normalcnt = co.cnt;
			double normaldepth = key.getDepth(normalcnt);
			double normaldepthAdj = normaldepth * normal.millionAdjust();
			long tumorcnt = co2.cnt;
			double tumordepth = key.getDepth(tumorcnt);
			double tumordepthAdj = tumordepth * tumor.millionAdjust();
			double tnratio = tumordepthAdj / normaldepthAdj;
			//
			if (!Double.isInfinite(normaldepthAdj)
					&& !Double.isNaN(normaldepthAdj)) {

				ss.addValue(normaldepthAdj);
			}
			if (!Double.isInfinite(tnratio) && !Double.isNaN(tnratio)) {
				tnss.addValue(tnratio);
			}
		}
		//

		ite = m1.keySet().iterator();
		int cnt = 0;
		String chrom = "";
		while (ite.hasNext()) {
			cnt++;
			CapInterval key = ite.next();
			String _chrom = key.chr;
			chrSet.add(_chrom);
			if (!chrom.equals(_chrom)) {
				if (!list.isEmpty()) {
					plist.add(list);
				}
				list = new ArrayList<WaveletIF>();
			}
			chrom = _chrom;
			Coverage.Container co = m1.get(key);
			Coverage.Container co2 = m2.get(key);
			long normalcnt = co.cnt;
			double normaldepth = key.getDepth(normalcnt);
			double normaldepthAdj = normaldepth * normal.millionAdjust();
			long tumorcnt = co2.cnt;
			double tumordepth = key.getDepth(tumorcnt);
			double tumordepthAdj = tumordepth * tumor.millionAdjust();
			double tnratio = tumordepthAdj / normaldepthAdj;
			boolean notSexChorm = (notSexChrom(chrom));
			tnratio = tnAdjust(tnratio, normalcnt, tumorcnt, normaldepthAdj,
					tumordepthAdj, ss, tnss, key, notSexChorm);

			System.out.println(normaldepthAdj + "\t" + tumordepthAdj);
			CNVInfo cnvinfo = new CNVInfo(normalcnt, normaldepth,
					normaldepthAdj, tumorcnt, tumordepth, tumordepthAdj,
					tnratio);
			key.setCNVInfo(cnvinfo);
			list.add(key);

		}
		plist.add(list);
		capinterval = plist;
		capintervalinit = true;
		return plist;

	}

	private boolean notSexChrom(String chrom) {

		chrom = chrom.toUpperCase();
		if (chrom.contains("X")) {
			return false;
		}
		if (chrom.contains("Y")) {
			return false;
		}
		return true;
	}

	private double tnAdjust(double tnratio, long normalcnt, long tumorcnt,
			double normaldepthAdj, double tumordepthAdj, SummaryStatistics ss,
			SummaryStatistics tnss, CapInterval key, boolean notSexChorm) {

		//double mean = tnss.getMean();
		if (lowCnt(normalcnt) || Double.isInfinite(tnratio)
				|| Double.isNaN(tnratio)) {
			
			
			//tnratio = tumordepthAdj / ss.getMean();
			tnratio = 1;
			//
//			if (Double.isInfinite(tnratio) || Double.isNaN(tnratio)) {
//				tnratio = mean;
//			}

		}
		//
		if (useAvearageNormal && notSexChorm) {

			//
			float nad = key.getNormalAveDepth();
			// nad = nad*100;
			if (nad > 0) {
				double fromnorm = tumordepthAdj / nad;
				if (fromnorm > 10 || fromnorm < 0.1) {

				} else {
					tnratio = fromnorm;
				}

			}

		}

		//

		// if (outof2sd(normaldepthAdj, ss)) {
		// double tn1 = tumordepthAdj / ss.getMean();
		// if (!Double.isInfinite(tn1) && !Double.isNaN(tn1)) {
		// if (Math.abs(mean - tn1) < Math.abs(mean - tnratio)) {
		// if (Math.abs(tnratio - tn1) > 0.1) {
		// return tn1;
		// }
		// }
		// }
		//
		// }
		//
		return tnratio;
	}

	private boolean outof2sd(double val, SummaryStatistics ss) {

		double diff = ss.getMean() - val;
		boolean loerend = (diff > (ss.getStandardDeviation() * 2));
		return loerend;
	}

	public List<String> getChromList() {
		List<String> chromList = new ArrayList<String>();
		for (String s : chrSet) {
			chromList.add(s);
		}
		return chromList;
	}

	private boolean lowCnt(long cnt) {

		return cnt < KarkinosProp.mincoverbase;
	}

	private boolean lowCnt2(long cnt) {

		return cnt < (KarkinosProp.mincoverbase * 10);
	}

	double baselineLOHEstimate;

	public double getBaselineLOHEstimate() {
		return baselineLOHEstimate;
	}

	public void setBaselineLOHEstimate(double _baselineLOHEstimate) {
		baselineLOHEstimate = _baselineLOHEstimate;
	}

	Map<Float, DataHolderByCN> map = null;
	AnalyseDist as = null;

	public AnalyseDist getAnalyseDist() {
		if (as == null) {
			//
			as = new AnalyseDist();
			as.analyseDist(this);
		}
		return as;
	}

	float fixtc = -1;

	public void setFixtc(float fixtc) {
		this.fixtc = fixtc;
	}

	double tumorratioFiitiingMatrix = 0;
	int baseploidy = 2;

	public void setBaseploidy(int _baseploidy) {
		baseploidy = _baseploidy;
	}

	public int getBaseploidy() {
		return baseploidy;
	}

	public float getTumorRatio() {

		if (fixtc > 0) {
			return fixtc;
		}

		float tratio = 1f;

		if (as != null) {

			tratio = as.getTumorratio();

		}
		if (tumorratioFiitiingMatrix > 0) {
			tratio = (float) tumorratioFiitiingMatrix;
			if (as != null) {
				if (as.getTcflg() == 3) {
					if (as.getTumorratio() > (tratio * 2)) {
						tratio = as.getTumorratio();
					}
				} else if (as.getTcflg() == 1) {

					if (as.getTumorratioFromLOH().getNumber() > 300) {
						tratio = as.getTumorratio();
					}
				}

			}
			boolean notenoughLOHSNP = as.getTumorratioFromLOH().getNumber() < 300;
			if ((tratio == 1) || (baseploidy == 4) || ((tratio < 0.4)&&notenoughLOHSNP)) {
				tratio = (float) tumorratioFiitiingMatrix;
			}
		} else {
			// no CNV peak
			if (as != null) {
				try {
					tratio = (float) as.getTumorratioFromSomatic()
							.getObservedratio();
				} catch (Exception ex) {

				}
				if (tratio > 0.2) {
					tratio = KarkinosProp.mintumorpurity;
				}
			}
		}
		if (tratio < KarkinosProp.mintumorpurity) {
			tratio = KarkinosProp.mintumorpurity;
		}
		return tratio;
	}

	public double getTumorratioFiitiingMatrix() {
		return tumorratioFiitiingMatrix;
	}

	public void setTumorratioFiitiingMatrix(double tumorratioFiitiingMatrix) {
		this.tumorratioFiitiingMatrix = tumorratioFiitiingMatrix;
	}

	FunctionRegression gcfr = null;

	public void setGcFunctionRegression(FunctionRegression fr) {
		gcfr = fr;

	}

	public FunctionRegression getFunctionRegression() {
		return gcfr;
	}

	List<CopyNumberInterval> cniVaridateList;

	public void setVaridatedCNlist(List<CopyNumberInterval> cniList) {
		cniVaridateList = cniList;
	}

	public List<CopyNumberInterval> getCniVaridateList() {
		return cniVaridateList;
	}



}
