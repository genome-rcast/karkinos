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
package jp.ac.utokyo.rcast.karkinos.filter;

import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.INFO_LOW_refOddsRatio;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.INFO_LOW_tumorOddsRatio;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;

public class FilterResult implements java.io.Serializable {

	// public static int SUPPORTED_BY_ONEDirection = 1;
	// public static int CONTAIN_INDEL_MISMATCH = 2;
	// public static int CONTAIN_MISMATCH = 3;
	// public static int READSENDSONLY = 4;
	public static final int PASS_FILTER = 0;
	public static final int SNP = 1;
	// move snp to INFO
	// public static final int SNP = 104;

	public static final int Low_tumor_adjustedRatio = 2;
	public static final int Low_tumor_adjustedReads = 3;
	public static final int High_normal_adjustedRatio = 4;
	public static final int LOW_adjusted_reads = 24;

	public static final int Low_mappability = 5;
	public static final int Low_complexty = 6;
	public static final int Low_SNPQual = 7;
	public static final int Low_MapQuality = 8;
	public static final int Basian_Filteling = 9;
	public static final int SUPPORTED_READSNG = 19;
	public static final int SUPPORTED_BY_ONEDirection = 11;
	public static final int CONTAIN_Reccurent_MISMATCH = 12;
	public static final int READSENDSONLY = 13;
	public static final int TooManyMismatchReads = 14;
	public static final int MutationAtSameCycle = 15;
	public static final int NEARINDEL = 16;
	public static final int HighNormalRatio = 17;
	public static final int noStrandSpecific = 18;
	//
	public static final int LOW_PROPER = 20;
	public static final int NOISE_IN_NORMAL = 25;

	public static final int highSecondMutatedAllele = 23;
	public static final int illuminaSpecific = 26;
	public static final int softClip = 27;
	public static final int TNQualityDiff = 28;

	public static final int FisherTestFail = 21;
	public static final int HighLowQualReads = 22;

	public static final int LOGn = 31;
	public static final int LOGt = 32;
	
	public static final int INFO_ffpe = 135;
	public static final int INFO_oxoG = 136;

	public static final int INFO_SUPPORTED_BY_ONEDirection = 101;
	public static final int INFO_Maskedrepeat_BY_RepeatMasker = 102;
	public static final int INFO_COSMIC = 103;
	public static final int INFO_COSMIC_Validate = 104;
	public static final int INFO_INIT_Filter_FAIL = 105;
	public static final int INFO_AllelicInfoAvailable = 106;
	public static final int INFO_SNP_flg = 107;

	public static final int INFO_LOW_refOddsRatio = 201;
	public static final int INFO_LOW_tumorOddsRatio = 202;
	public static final int INFO_adjustAlleleFreq = 203;
	public static final int INFO_adjustLowdepth = 204;
	public static final int INFO_minimumSupportReads = 205;

	// bw.append("##FILTER=<ID=qf,Description=\"Quality below threshold\">");
	// bw.append("##FILTER=<ID=bf,Description=\"Bayesian filterling\">");
	// bw.append("##FILTER=<ID=snp,Description=\"dbSNP snp\">");
	// bw.append("##FILTER=<ID=ma,Description=\"Low mappability \">");
	// bw.append("##FILTER=<ID=entropy,Description=\"Low complexty\">");
	// bw.append("##FILTER=<ID=srd,Description=\"Support reads have only one direction\">");
	// bw.append("##FILTER=<ID=srm,Description=\"Only support reads have another mismatch\">");
	// bw.append("##FILTER=<ID=sre,Description=\"mutation at reads ends\">");

	String before10;
	String after10;

	boolean typicalSysErr;

	public boolean isTypicalSysErr() {
		return typicalSysErr;
	}

	public void setTypicalSysErr(boolean typicalSysErr) {
		this.typicalSysErr = typicalSysErr;
	}

	public String getBefore10() {
		return before10;
	}

	public void setBefore10(String before10) {
		this.before10 = before10;
	}

	public String getAfter10() {
		return after10;
	}

	public void setAfter10(String after10) {
		this.after10 = after10;
	}

	public String getFilterStr(Set<Integer> filter) {

		StringBuffer sb = new StringBuffer();
		for (int n : filter) {

			if (n > 100)
				continue;
			if (sb.length() > 0) {
				sb.append(",");
			}
			sb.append(_getFilterStr(n));
		}
		return sb.toString();

	}

	// public static final int SNP = 1;
	// public static final int Low_mappability = 5;
	// public static final int Low_complexty = 6;
	// public static final int Low_SNPQual = 7;
	// public static final int Basian_Filteling = 9;
	// public static final int SUPPORTED_READSNG = 19;
	// public static final int SUPPORTED_BY_ONEDirection = 11;
	// public static final int CONTAIN_Reccurent_MISMATCH = 12;
	// public static final int READSENDSONLY = 13;
	// public static final int TooManyMismatchReads = 14;
	// public static final int MutationAtSameCycle = 15;
	public String _getFilterStr(int flg) {

		//
		switch (flg) {
		case SNP:
			return "snp";
		case Basian_Filteling:
			return "bf";
		case Low_SNPQual:
			return "qf";
		case Low_mappability:
			return "ma";
		case Low_MapQuality:
			return "mq";
		case Low_complexty:
			return "entropy";
		case READSENDSONLY:
			return "sre";
		case CONTAIN_Reccurent_MISMATCH:
			return "srm";
		case SUPPORTED_BY_ONEDirection:
			return "srd";
		case TooManyMismatchReads:
			return "mmt";
		case MutationAtSameCycle:
			return "scc";
		case FisherTestFail:
			return "fisher";
		case LOGn:
			return "logn";
		case LOGt:
			return "logt";
		case NEARINDEL:
			return "near_indel";
		case Low_tumor_adjustedRatio:
			return "low_adj_ratio";
		case Low_tumor_adjustedReads:
			return "low_adj_reads";
		case HighLowQualReads:
			return "high_lqr";
		case High_normal_adjustedRatio:
			return "high_TCn";
		case LOW_adjusted_reads:
			return "lowTCreads";
		case INFO_adjustLowdepth:
			return "lowTCreads2";
		case SUPPORTED_READSNG:
			return "srng";
		case HighNormalRatio:
			return "highNratio";
		case noStrandSpecific:
			return "nss";
		case highSecondMutatedAllele:
			return "hsa";
		case illuminaSpecific:
			return "illuminaSys";
		case softClip:
			return "softClip";
		case TNQualityDiff:
			return "tnQdiff";
			// public static final int LOW_PROPER = 20;
			// public static final int NOISE_IN_NORMAL = 21;
		case LOW_PROPER:
			return "lowProper";
		case NOISE_IN_NORMAL:
			return "noise_in_normal";
		
			// High_normal_adjustedRatio
			// Low_tumor_adjustedRatio
			// Low_tumor_adjustedReads
			// LOW_adjusted_reads

			// case High_normal_adjustedRatio:
			// return "high_adj_ratio";
		}

		if (flg != 0) {
			return String.valueOf(flg);
		}

		return "PASS";
	}

	Set<Integer> filter = null;

	public Set<Integer> getPassFilterFlg() {
		if (filter == null) {
			filter = _getPassFilterFlg();
		}
		return filter;
	}

	public Set<Integer> getPassFilterFlgForce() {
		filter = _getPassFilterFlg();
		return filter;
	}

	public Set<Integer> _getPassFilterFlg() {

		Set<Integer> filter = new LinkedHashSet<Integer>();
		if (supportReadsFlgs == null) {
			supportReadsFlgs = new HashSet<Integer>();
		}
		boolean cosmicvalid = false;
		//
		if (normalVrcnt > 0 && tumorAF > 0) {
			//
			// double ratio = (double)((double)tumorVrcnt/(double)normalVrcnt);
			double ratio = normalAF / tumorAF;
			if (ratio > KarkinosProp.minTumorNormalRatio) {
				filter.add(HighNormalRatio);
			}

		}
		boolean lowsupport = (numSupportRead <= 6);

		if (dbSNPbean != null) {

			boolean dbSNP = dbSNPbean.getMode() == DbSNPAnnotation.MODEdbSNP;
			boolean onekg = (dbSNPbean.getMode() == DbSNPAnnotation.MODE1000g);
			boolean exonSNP = (dbSNPbean.getMode() == DbSNPAnnotation.MODEexonSNP);

			boolean lowdapthLow = normaldepth < KarkinosProp.low_normal_depth_thresLow;
			boolean lowdapthHigh = normaldepth < KarkinosProp.low_normal_depth_thresHigh;

			// if (dbSNPbean.isCosmic()) {
			//
			// // filter.add(INFO_COSMIC);
			// supportReadsFlgs.add(INFO_COSMIC);
			// if (dbSNPbean.isValid()) {
			// supportReadsFlgs.add(INFO_COSMIC_Validate);
			// cosmicvalid = dbSNPbean.isCosmicHigh();
			// }
			//
			// if(dbSNP&&highisoform){
			// filter.add(SNP);
			// }
			//
			// } else {
			//
			//
			// boolean onekg = (dbSNPbean.getMode() ==
			// DbSNPAnnotation.MODE1000g);
			// boolean exonSNP = (dbSNPbean.getMode() ==
			// DbSNPAnnotation.MODEexonSNP);
			// boolean lowdapth = normaldepth <
			// KarkinosProp.low_normal_depth_thres;
			//
			// if (dbSNPbean.isValid()) {
			//
			// filter.add(SNP);
			//
			// } else if (lowdapth || highisoform) {
			//
			// // dbSNP and low normal
			// if (dbSNP || onekg || exonSNP) {
			//
			// filter.add(SNP);
			// }
			// }
			// }
			// policy change 20130308

			if (dbSNPbean.isCosmic()) {

				// filter.add(INFO_COSMIC);
				supportReadsFlgs.add(INFO_COSMIC);
				if (dbSNPbean.isCosmicvalid()) {
					supportReadsFlgs.add(INFO_COSMIC_Validate);
					// cosmicvalid = dbSNPbean.isCosmicHigh();
					cosmicvalid = true;
				}

			} else if (dbSNP || onekg || exonSNP) {

				if (indel) {
					filter.add(SNP);
				}
				boolean validated = dbSNPbean.isValid();
				// if (lowdapthHigh) {
				// if (validated) {
				// filter.add(SNP);
				// }
				// }
				// change v4.1.10
				if (validated) {
					filter.add(SNP);
				}

				if (lowsupport || lowdapthLow || uncertainDelation
						|| (lowsupport && (tumorAF <= 0.2))) {
					filter.add(SNP);
				}

			}

			if (!dbSNPbean.isCosmic()) {
				if (dbSNPbean.isValid()) {

					supportReadsFlgs.add(INFO_SNP_flg);

				}
			}

		}
		//

		if (numSupportRead <= KarkinosProp.minsupportreads) {

			supportReadsFlgs.add(INFO_minimumSupportReads);
			if (adjustedTumorAllereFreq < 0.3 || normalVrcnt > 0) {
				filter.add(Low_tumor_adjustedRatio);
			}
			boolean lowdepth = normaldepth < KarkinosProp.low_normal_depth_thresLow;
			if (lowdepth) {
				filter.add(LOW_adjusted_reads);
			}
		}

		if (initFilterFalse) {

			supportReadsFlgs.add(INFO_INIT_Filter_FAIL);

		}

		// take out side of !cosmicvalid block 20130308
		if (highnormalAfterTCadjust) {
			filter.add(High_normal_adjustedRatio);
		}

		if (!cosmicvalid) {

			// except indel case
			// if (!indel) {
			// if (adjustedTumorAllereFreq <
			// KarkinosProp.mintumorratioForFilter1) {
			// filter.add(Low_tumor_adjustedRatio);
			// }
			// }

			double mintratio = KarkinosProp.mintumorratio;
			if (initFilterFalse) {
				mintratio = KarkinosProp.mintumorratioForResqued;
			}
			if (indel) {
				mintratio = (mintratio * 0.5);
			}
			if (adjustedTumorAllereFreq < mintratio) {
				// except indel case
				supportReadsFlgs.add(INFO_adjustAlleleFreq);
			}

			if (lowtumorAfterTCadjust) {

				filter.add(Low_tumor_adjustedReads);

			}

			if (lowdepthafterTCadjust) {

				filter.add(LOW_adjusted_reads);

			}
			if (lowdepthafterTCadjustInfo) {

				supportReadsFlgs.add(INFO_adjustLowdepth);

			}
		}

		if (!fisherTestSignif) {
			filter.add(FisherTestFail);
		}

		if (mapQ < KarkinosProp.minMapQ) {
			filter.add(Low_MapQuality);
		}
		// /

		// v 4.1.10 comment entropy depth
		// if (tumordepth < KarkinosProp.entropyDepth) {
		if (ignoreEntroptCheck == false && (tumorAF <= 0.4)) {
			if (seqEntropy < KarkinosProp.minEntropy) {
				filter.add(Low_complexty);
			}
		}

		// }
		if (lowsupport) {
			if (seqEntropy < 1.0) {
				filter.add(Low_complexty);
			}
		}

		if (lowsupport) {

			if (typicalSysErr) {
				if (phredbq < 120) {
					filter.add(illuminaSpecific);
				}
			}

		}

		// filter 2013.10.02
		if (mappability < KarkinosProp.minMappability) {
			filter.add(Low_mappability);
		}

		if (!indel) {

			if (phredbq < KarkinosProp.minPhredQual) {
				filter.add(Low_SNPQual);
			}

			// if (bresult != null) {
			// if (bresult.isResult()) {
			// filter.add(Basian_Filteling);
			// }
			// }
			// basian filters
			boolean lowthresn = lognAjusted < KarkinosProp.LognThres;
			boolean lowdepthn = normaldepth < KarkinosProp.baysianFilterdepth;
			if (lowthresn) {
				if (lowdepthn || highisoform) {

					supportReadsFlgs.add(INFO_LOW_refOddsRatio);
				}
			}

			// boolean lowthrest = logtAjusted < KarkinosProp.LogtThres;
			// boolean logtlow = logt <= 2;
			// boolean lowdeptht = tumordepth < KarkinosProp.baysianFilterdepth;
			// boolean lowratio = adjustedTumorAllereFreq < 0.3;
			// if (lowthrest || logtlow) {
			//
			// if ((lowdeptht && lowratio) || highisoform) {
			// supportReadsFlgs.add(INFO_LOW_tumorOddsRatio);
			// }
			//
			// }
			// check second mutated allele
			double sAFratio = secondAF / tumorAF;
			if (sAFratio > KarkinosProp.minSecondMutationRelativeRatio) {
				filter.add(highSecondMutatedAllele);
			}

		}

		for (int flg : supportReadsFlgs) {
			if (flg < 100) {
				
				if(flg==CONTAIN_Reccurent_MISMATCH){
					if(cosmicvalid){
						continue;
					}
				}				
				filter.add(flg);
			}
		}

		if (filter.isEmpty()) {
			filter.add(PASS_FILTER);
		}

		return filter;
	}

	BfilterResult bresult;

	public BfilterResult getBresult() {
		return bresult;
	}

	public void setBresult(BfilterResult bresult) {
		this.bresult = bresult;
	}

	DbSNPBean dbSNPbean;

	public DbSNPBean getDbSNPbean() {
		return dbSNPbean;
	}

	public void setDbSNPbean(DbSNPBean dbSNPbean) {
		this.dbSNPbean = dbSNPbean;
	}

	double pvalFisher;
	boolean fisherTestSignif;
	double logt = 0;
	double logn = 0;
	double lognAjusted = 0;
	boolean uncertainDelation = false;
	float ntQualityDiff = 0;

	public void setNtQualityDiff(float ntQualityDiff) {
		this.ntQualityDiff = ntQualityDiff;
	}

	public void setUncertainDelation(boolean uncertainDelation) {
		this.uncertainDelation = uncertainDelation;
	}

	public void setHighisoform(boolean highisoform) {
		this.highisoform = highisoform;
	}

	boolean ignoreEntroptCheck;
	boolean highnormalAfterTCadjust;
	boolean lowtumorAfterTCadjust;
	boolean lowdepthafterTCadjust;
	int numSupportRead;
	boolean highisoform = false;
	int normalVrcnt = 0;
	int tumorVrcnt = 0;
	double normalAF = 0;
	double tumorAF = 0;
	float supportreadsBAlleleFeqquency = 0;
	float refreadsBAlleleFeqquency = 0;
	char secondAllele;
	float secondAF;

	public char getSecondAllele() {
		return secondAllele;
	}

	public void setSecondAllele(char secondAllele) {
		this.secondAllele = secondAllele;
	}

	public float getSecondAF() {
		return secondAF;
	}

	public void setSecondAF(double d) {
		this.secondAF = (float) d;

	}

	public float getSupportreadsBAlleleFeqquency() {
		return supportreadsBAlleleFeqquency;
	}

	public void setSupportreadsBAlleleFeqquency(
			float supportreadsBAlleleFeqquency) {
		this.supportreadsBAlleleFeqquency = supportreadsBAlleleFeqquency;
	}

	public float getRefreadsBAlleleFeqquency() {
		return refreadsBAlleleFeqquency;
	}

	public void setRefreadsBAlleleFeqquency(float refreadsBAlleleFeqquency) {
		this.refreadsBAlleleFeqquency = refreadsBAlleleFeqquency;
	}

	public void setNormalAF(double normalAF) {
		this.normalAF = normalAF;
	}

	public void setTumorAF(double tumorAF) {
		this.tumorAF = tumorAF;
	}

	public void setNormalVrcnt(int normalVrcnt) {
		this.normalVrcnt = normalVrcnt;
	}

	public void setTumorVrcnt(int tumorVrcnt) {
		this.tumorVrcnt = tumorVrcnt;
	}

	public int getNumSupportRead() {
		return numSupportRead;
	}

	public void setNumSupportRead(int numSupportRead) {
		this.numSupportRead = numSupportRead;
	}

	boolean lowdepthafterTCadjustInfo;

	public boolean isLowdepthafterTCadjustInfo() {
		return lowdepthafterTCadjustInfo;
	}

	public void setLowdepthafterTCadjustInfo(boolean lowdepthafterTCadjustInfo) {
		this.lowdepthafterTCadjustInfo = lowdepthafterTCadjustInfo;
	}

	public boolean isLowdepthafterTCadjust() {
		return lowdepthafterTCadjust;
	}

	public void setLowdepthafterTCadjust(boolean lowdepthafterTCadjust) {
		this.lowdepthafterTCadjust = lowdepthafterTCadjust;
	}

	public boolean isLowtumorAfterTCadjust() {
		return lowtumorAfterTCadjust;
	}

	public void setLowtumorAfterTCadjust(boolean lowtumorAfterTCadjust) {
		this.lowtumorAfterTCadjust = lowtumorAfterTCadjust;
	}

	public boolean isHighnormalAfterTCadjust() {
		return highnormalAfterTCadjust;
	}

	public void setHighnormalAfterTCadjust(boolean highnormalAfterTCadjust) {
		this.highnormalAfterTCadjust = highnormalAfterTCadjust;
	}

	double adjustedTumorAllereFreq;
	int adjustedtumortotalreads;

	public double getAdjustedTumorAllereFreq() {
		return adjustedTumorAllereFreq;
	}

	public void setAdjustedTumorAllereFreq(double adjustedTumorAllereFreq) {
		this.adjustedTumorAllereFreq = adjustedTumorAllereFreq;
	}

	public int getAdjustedtumortotalreads() {
		return adjustedtumortotalreads;
	}

	public void setAdjustedtumortotalreads(int adjustedtumortotalreads) {
		this.adjustedtumortotalreads = adjustedtumortotalreads;
	}

	public void setIgnoreEntroptCheck(boolean ignoreEntroptCheck) {
		this.ignoreEntroptCheck = ignoreEntroptCheck;
	}

	public double getLognAjusted() {
		return lognAjusted;
	}

	public void setLognAjusted(double lognAjusted) {
		this.lognAjusted = lognAjusted;

	}

	double logtAjusted = 0;

	public double getLogtAjusted() {
		return logtAjusted;
	}

	public void setLogtAjusted(double adjustedLogt) {

		this.logtAjusted = adjustedLogt;

	}

	public double getLogt() {
		return logt;
	}

	public void setLogt(double logt) {
		this.logt = logt;
	}

	public double getLogn() {
		return logn;
	}

	public void setLogn(double logn) {
		this.logn = logn;
	}

	public double getPvalFisher() {
		return pvalFisher;
	}

	public void setPvalFisher(double pvalFisher) {
		this.pvalFisher = pvalFisher;
	}

	float pval4FiserDirectional = 0f;

	public float getPval4FiserDirectional() {
		return pval4FiserDirectional;
	}

	public void setPval4FiserDirectional(float pval4FiserDirectional) {
		this.pval4FiserDirectional = pval4FiserDirectional;
	}

	public boolean isFisherTestSignif() {
		return fisherTestSignif;
	}

	public void setFisherTestSignif(boolean fisherTestSignif) {
		this.fisherTestSignif = fisherTestSignif;
	}

	double seqEntropy;
	boolean repeat;
	float mappability;
	Set<Integer> supportReadsFlgs;
	float phredbq;
	boolean indel;
	float mapQ;
	int normaldepth;
	int tumordepth;

	public int getNormaldepth() {
		return normaldepth;
	}

	public void setNormaldepth(int normaldepth) {
		this.normaldepth = normaldepth;
	}

	public int getTumordepth() {
		return tumordepth;
	}

	public void setTumordepth(int tumordepth) {
		this.tumordepth = tumordepth;
	}

	public float getMapQ() {
		return mapQ;
	}

	public void setMapQ(float mapQ) {
		this.mapQ = mapQ;
	}

	public boolean isIndel() {
		return indel;
	}

	public void setIndel(boolean indel) {
		this.indel = indel;
	}

	public float getPhredbq() {
		return phredbq;
	}

	public void setPhredbq(float phredbq) {
		this.phredbq = phredbq;
	}

	public double getSeqEntropy() {
		return seqEntropy;
	}

	public void setSeqEntropy(double seqEntropy) {
		this.seqEntropy = seqEntropy;
	}

	public boolean isRepeat() {
		return repeat;
	}

	public void setRepeat(boolean reapeat) {
		this.repeat = reapeat;
	}

	public float getMappability() {
		return mappability;
	}

	public void setMappability(float mappability) {
		this.mappability = mappability;
	}

	public void setSupportReadsFlgs(Set<Integer> supportReadsFlgs) {
		this.supportReadsFlgs = supportReadsFlgs;
	}

	public boolean isPassFilter() {

		Set<Integer> filter = getPassFilterFlg();
		if (filter == null || filter.isEmpty()) {
			return true;
		}
		// 2013/01/18
		// if (dbSNPbean != null) {
		// if (dbSNPbean.isCosmic()) {
		//
		// return true;
		// }
		// }
		return filter.contains(PASS_FILTER);

	}

	public Set<Integer> getInfoFlg() {
		Set<Integer> info = new HashSet<Integer>();
		for (int flg : supportReadsFlgs) {
			if (flg >= 100) {
				info.add(flg);
			}
		}
		return info;
	}

	boolean initFilterFalse = false;

	public void setInitFilterFalse(boolean initFilterFalse) {

		this.initFilterFalse = initFilterFalse;

	}

}
