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

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

public class PileUPResult implements java.io.Serializable {

	public void setRatiocalclated(boolean ratiocalclated) {
		this.ratiocalclated = ratiocalclated;
	}

	boolean diff = false;
	int[] seqCounter = new int[5];
	double[] expected = new double[5];
	double[] likehood = new double[5];
	double[][] likehoodp2x = new double[5][5];
	double[] phredQual = new double[5];
	double[] mapqualRMS = new double[5];
	double[] phredQualRMS = new double[5];
	double refp = 0;
	double mutatep = 0;

	IndelInfo indelinfo = null;
	Map<String, Counter> insersionmap = null;
	Map<Integer, Counter> delmap = null;
	boolean indel = false;
	int totalcnt = 0;
	double mapqualRMS4Indel = 0;
	float lowqualratio = 0f;

	public void clear() {

		indelinfo = null;
		diff = false;
		ratiocalclated = false;
		maxidx = 0;
		genomeR = 0;
		total = 0;
		insersionmap = null;
		delmap = null;
		indel = false;
		indelcnt = 0;
		noindelcnt = 0;
		insersion = false;
		totalcnt = 0;
		mapqualRMS4Indel = 0;
		double refp = 0;
		double mutatep = 0;
		lowqualratio = 0f;

		for (int n = 0; n < 5; n++) {
			seqCounter[n] = 0;
			likehood[n] = 0;
			phredQual[n] = 0f;
			phredQualRMS[n] = 0;
			expected[n] = 0;
			mapqualRMS[n] = 0;
			for (int m = 0; m < 5; m++) {
				likehoodp2x[n][m] = 0;
			}
		}

	}

	public float getLowqualratio() {
		return lowqualratio;
	}

	public void setLowqualratio(float lowqualratio) {
		this.lowqualratio = lowqualratio;
	}

	boolean insersion = false;

	public String getIndelStr() {

		if (indel) {
			if (isInsersion()) {
				StringBuffer sb = new StringBuffer();
				Set<String> s = insersionmap.keySet();
				for (String ss : s) {
					sb.append(ss + "\t");
				}
				return sb.toString();
			} else {
				StringBuffer sb = new StringBuffer();
				Set<Integer> i = delmap.keySet();
				for (Integer in : i) {
					sb.append(+in + "\t");
				}
				return sb.toString();
			}
		}
		return "";
	}

	public double getTotal() {
		return totalcnt;
	}

	public char getGenomeR() {
		return genomeR;
	}

	private int getMaxidx() {

		double maxdd = 0;
		int idx = 0;
		int maxidx = 0;
		int gidx = seqALL.indexOf(genomeR);

		for (double d : expected) {

			if (gidx == idx) {
				idx++;
				continue;
			}
			if (d > maxdd) {
				maxdd = d;
				maxidx = idx;
			}
			idx++;
		}
		return maxidx;
	}

	private int getSecondMutationidx(){
		
		int gidx = seqALL.indexOf(genomeR);
		int maxidx = getMaxidx();
		if(indel){
			return maxidx;
		}
		
		int sndidx = 0;
		double maxdd = 0;
		int idx = 0;

		for (double d : expected) {

			if ((gidx == idx) || (maxidx == idx)) {
				idx++;
				continue;
			}
			if (d > maxdd) {
				maxdd = d;
				sndidx = idx;
			}
			idx++;
		}
		return sndidx;
		
	}
	
	public void setIndel(boolean indel) {
		this.indel = indel;
		ratiocalclated = false;
	}

	
	public int getRefCnt() {
		
		int gidx = seqALL.indexOf(genomeR);
		return seqCounter[gidx];
		
	}

	public int getAltCnt() {
		if (indel) {
			return indelcnt();
		}else{
			return seqCounter[getMaxidx()];
		}
	}

	public int indelcnt() {
		int indelcnt = 0;
		if (delmap != null) {
			indelcnt = getMax(delmap);
		}
		if (insersionmap != null) {
			int indelcnt2 = getMax(insersionmap);
			if(indelcnt2>indelcnt)indelcnt = indelcnt2;
		}
		return indelcnt;
	}

	public char getALT() {

		return seqALL.charAt(getMaxidx());

	}
	
	public char getSecondALT(){
		
		return seqALL.charAt(getSecondMutationidx());
		
	}

	public float getPhred() {

		return (float) phredQual[getMaxidx()] * 10;

	}

	public float getBQrms() {
		if (indel) {
			return 0f;
		}
		int mi = getMaxidx();
		int count = seqCounter[mi];
		double pq = phredQualRMS[mi];
		double d = pq / (double) count;
		return (float) Math.sqrt(d);

	}

	public float getMQrms() {

		if (indel) {
			double d = mapqualRMS4Indel / (double) indelcnt;
			return (float) Math.sqrt(d);
		}
		int mi = getMaxidx();
		int count = seqCounter[mi];
		double pq = mapqualRMS[mi];
		double d = pq / (double) count;
		return (float) Math.sqrt(d);
	}

	private static double getBinomial(double p_hat, int cnt, int total)
			throws MathException {

		BinomialDistribution bd = new BinomialDistributionImpl(total, p_hat);
		double p = bd.cumulativeProbability(cnt);
		return p;
	}

	public class Counter implements java.io.Serializable {
		int n = 0;

		void inc() {
			n++;
		}
	}

	public static final String seqALL = "ATCGN";
	int indelcnt = 0;
	private int noindelcnt = 0;

	public void setBaseAndQual(char ch, byte qual, int mapq, IndelInfo indelinfo) {

		totalcnt++;
		if (indelinfo.indel) {
			indelcnt++;
			indel = true;
			boolean del = (indelinfo.insersion == null);
			if (del) {
				int len = indelinfo.length;
				if (delmap == null)
					delmap = new HashMap<Integer, Counter>();

				Counter ct = null;
				if (delmap.containsKey(len)) {
					ct = delmap.get(len);
				} else {
					ct = new Counter();
					delmap.put(len, ct);
				}
				ct.inc();
			} else {
				if (insersionmap == null)
					insersionmap = new HashMap<String, Counter>();
				Counter ct = null;
				if (insersionmap.containsKey(indelinfo.insersion)) {
					ct = insersionmap.get(indelinfo.insersion);
				} else {
					ct = new Counter();
					insersionmap.put(indelinfo.insersion, ct);
				}
				insersion = true;
				ct.inc();
			}
			mapqualRMS4Indel = mapqualRMS4Indel + Math.pow(mapq, 2);
		} else {
			noindelcnt++;
		}

		if (ch > 0) {
			//
			int idx = seqALL.indexOf(ch);
			if (idx < 0) {
				// something wrong
			} else {
				seqCounter[idx] = seqCounter[idx] + 1;
				double qual0 = (int) qual & 0xFF;

				phredQualRMS[idx] = phredQualRMS[idx] + Math.pow(qual0, 2);
				mapqualRMS[idx] = mapqualRMS[idx] + Math.pow(mapq, 2);

				qual0 = qual0 * 0.1;
				phredQual[idx] = phredQual[idx] + qual0;
				double pNomatch = (1 / Math.pow(10, qual0));
				double pmathch = 1 - pNomatch;

				for (int n = 0; n < 4; n++) {

					if (n == idx) {
						expected[n] = expected[n] + pmathch;
						likehoodp2x[idx][n] = likehoodp2x[idx][n] + pmathch;
						// double logpmatch = Math.log10(pmathch);
						// likehood[n] = likehood[n] + logpmatch;

					} else {
						expected[n] = expected[n] + (pNomatch / (double) 3);
						likehoodp2x[idx][n] = likehoodp2x[idx][n]
								+ (pNomatch / (double) 3);
						// double logpunmatch =Math.log10(pNomatch/(double)3);
						// likehood[n] = likehood[n] +logpunmatch;

					}

				}

			}

		}

	}

	public int getTotalcnt() {
		return totalcnt;
	}

	public void setTotalcnt(int totalcnt) {
		this.totalcnt = totalcnt;
	}

	public boolean isInsersion() {
		
		return getMax(insersionmap) > getMax(delmap);
		
	}

	public boolean isIndel() {
		return indel;
	}

	public Map<String, Counter> getInsersionmap() {
		return insersionmap;
	}

	public Map<Integer, Counter> getDelmap() {
		return delmap;
	}

	public boolean isDiff() {
		return diff;
	}

	public int[] getSeqCounter() {
		return seqCounter;
	}

	public double[] getPhredQual() {
		return phredQual;
	}

	public void setDiff(boolean _diff) {
		diff = _diff;
	}

	boolean ratiocalclated = false;
	float ratio = 0f;
	int maxidx = 0;
	double total = 0;

	private void calcRatio() {

		if (indel) {

			double indelr = calcRatio_Indel();
			double snvr = calcRatio_SNV();
//			if (snvr > indelr) {
//				ratio = (float) snvr;
//				indel = false;
//			} else {
				ratio = (float) indelr;
//			}
		} else {
			calcRatio_SNV();
		}
	}

	private double calcRatio_Indel() {

		int indelcnt = indelcnt();		
		ratio = (float) ((double) indelcnt / (double) (indelcnt + noindelcnt));
		return ratio;
	}

	public float getRefLogLikeHood() {

		int mi = getMaxidx();
		int refidx = seqALL.indexOf(genomeR);
		if ((refidx < 0) || (mi < 0)) {
			return 100f;
		}
		double refd = expected[refidx];
		double mid = expected[mi];
		double d = refd / mid;
		return (float) (Math.log10(d));

	}

	public float getMutationLogLikeHood() {

		int mi = getMaxidx();
		int refidx = seqALL.indexOf(genomeR);
		if ((refidx < 0) || (mi < 0)) {
			return 100f;
		}
		double refd = expected[refidx];
		double mid = expected[mi];
		double d = mid / refd;
		return (float) (Math.log10(d));

	}

	public float getMutateLogLikeHoodAmongMutation() {

		int mi = getMaxidx();
		int refidx = seqALL.indexOf(genomeR);
		if ((refidx < 0) || (mi < 0)) {
			return 100f;
		}
		double refd = likehoodp2x[mi][refidx];
		double mid = likehoodp2x[mi][mi];
		double d = mid / refd;
		return (float) (Math.log10(d));

	}

	public int[] getRefAltCnt() {

		int[] ret = new int[2];
		if (indel) {
			if (noindelcnt == 0 && indelcnt == 0) {
				noindelcnt = totalcnt;
			}
			ret[0] = noindelcnt;
			ret[1] = indelcnt;
		} else {
			int mi = getMaxidx();
			int refidx = seqALL.indexOf(genomeR);
			if (refidx == -1) {
				ret[0] = totalcnt;
				ret[1] = 0;
				return ret;
			}
			int refcnt = seqCounter[refidx];
			int count = seqCounter[mi];
			if (refcnt == 0 && count == 0) {
				refcnt = totalcnt;
			}
			ret[0] = refcnt;
			ret[1] = count;
		}
		return ret;

	}

	private double calcRatio_SNV() {

		int idx = 0;
		int idxref = seqALL.indexOf(genomeR);
		double max = 0;

		for (double d : expected) {

			if (idx != idxref) {
				if (max == 0)
					max = d;
				if (max < d) {
					max = d;
					maxidx = idx;
				}
			}
			total = total + d;
			idx++;
		}
		if (total != 0) {
			ratio = (float) ((double) max / (double) total);
		} else {
			ratio = 0;
		}
		return ratio;
	}

	public float getRatio() {
		if (!ratiocalclated || ratio == 0f) {
			calcRatio();
			ratiocalclated = true;
		}
		return ratio;
	}

	public double getSecondRatio() {

		int idx = 0;
		int idxsecond = getSecondMutationidx();
		double second = 0;
		double l_total = 0;
		double l_ratio = 0;
		
		for (double d : expected) {

			if (idx == idxsecond) {
				second = d;
			}
			l_total = l_total+ d;
			idx++;
		}
		if (l_total != 0) {
			l_ratio = (float) ((double) second / (double) l_total);
		} else {
			l_ratio = 0;
		}
		return l_ratio;
	}
	
	public boolean haveLowerRatio(float _ratio) {
		double ratio = getRatio();
		return ratio <= _ratio;
	}

//	public boolean haveLowerRatio(float _ratio, char alt) {
//		double ratio = calcRatio_SNV(alt);
//		return ratio <= _ratio;
//	}
//
//	private double calcRatio_SNV(char alt) {
//		try {
//			int idxref = seqALL.indexOf(genomeR);
//			int idxalt = seqALL.indexOf(alt);
//			if ((idxref < 0) || (idxalt < 0)) {
//				return getRatio();
//			}
//			double totall = 0;
//			double ratiol = 0;
//			for (double d : expected) {
//				totall  = totall + d;
//			}
//
//			if (totall!= 0) {
//				ratiol = (float) ((double) expected[idxalt] / (double) total);
//			} else {
//				ratiol = 0;
//			}
//			return ratiol;
//		} catch (Exception ex) {
//			return getRatio();
//		}
//	}

	public boolean haveHigherRatio(float _ratio) {
		double ratio = getRatio();
		return ratio >= _ratio;
	}

	char genomeR = 0;

	public void setGenomeRef(char _genomeR) {
		genomeR = Character.toUpperCase(_genomeR);
	}

	public PileUPResult getIndelCopy() {
		
		PileUPResult ret = new PileUPResult();
		ret.totalcnt = totalcnt;
		ret.indel = indel;
		ret.ratiocalclated = false;
		ret.noindelcnt = noindelcnt;
		ret.indelcnt = indelcnt;
		ret.mapqualRMS4Indel = mapqualRMS4Indel;

		if (delmap != null) {
			ret.delmap = new HashMap<Integer, Counter>();
			deepcopy(delmap, ret.delmap);
		}
		if (insersionmap != null) {
			ret.insersionmap = new HashMap<String, Counter>();
			deepcopy(insersionmap, ret.insersionmap);
		}

		return ret;
	}

	private void deepcopy(Map source, Map dist) {
		Set<Entry> es = source.entrySet();
		for (Entry e : es) {
			dist.put(e.getKey(), e.getValue());
		}

	}

	private int getMax(Map m) {
		if(m==null)return 0;
		Set<Entry> es = m.entrySet();
		int max = 0;
		for (Entry e : es) {
			int cnt = ((Counter) e.getValue()).n;
			if (cnt > max) {
				max = cnt;
			}
		}
		return max;
	}

	public String getInfoStr() {
		StringBuffer sb = new StringBuffer();
		for (double d : expected) {
			sb.append(d + "\t");
		}
		return sb.toString();
	}

}
