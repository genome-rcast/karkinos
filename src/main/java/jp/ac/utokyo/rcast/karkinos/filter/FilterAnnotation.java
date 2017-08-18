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

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.MappabilityAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.stats.Fisher;
import jp.ac.utokyo.rcast.karkinos.distribution.AnalyseDist;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.readssummary.SNPDepthCounter;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;
import jp.ac.utokyo.rcast.karkinos.utils.GenotypeKeyUtils;
import jp.ac.utokyo.rcast.karkinos.utils.NormalSNPCounter;
import jp.ac.utokyo.rcast.karkinos.utils.SNVHighDepthCounter;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class FilterAnnotation {

	String mappability = null;
	TwoBitGenomeReader tgr = null;
	MappabilityAnnotation mappabilityanno = null;
	String normalbam = null;
	String tumorbam = null;
	SupportReadsCheck srCheck = null;
	boolean checkMappability = false;
	GeneExons ge = null;
	NoiseAnalysis na = null;

	public FilterAnnotation(String mappability, TwoBitGenomeReader tgr, String normalbamf, String tumorbamf,
			DbSNPAnnotation dbAnno, GeneExons ge) throws IOException {
		this.mappability = mappability;
		mappabilityanno = new MappabilityAnnotation(mappability);
		checkMappability = mappabilityanno.isInit();
		this.tgr = tgr;
		this.normalbam = normalbamf;
		this.tumorbam = tumorbamf;
		srCheck = new SupportReadsCheck(normalbam, tumorbam, tgr, dbAnno);
		this.ge = ge;

	}

	public static int debugpos = 189868614;

	// 20130403 first filter just to use to find peak pulidity
	public void filterAnnotation1(DataSet dataset, ReadsSummary readsSummary) throws IOException {

		NormalSNPCounter nsnpC = new NormalSNPCounter(dataset);

		Map<String, Integer> snppos = new HashMap<String, Integer>();

		for (SNVHolder snv : dataset.getSnvlist()) {

			int flg = snv.getFlg();
			boolean check = (flg == PileUP.SomaticMutation || snv.getFlg() == PileUP.TumorINDEL);
			if (flg == PileUP.REGBOTH || flg == PileUP.NormalSNP || flg == PileUP.BothINDEL
					|| flg == PileUP.NormalINDEL) {
				String snpposS = snv.getChr() + "-" + snv.getPos();
				snppos.put(snpposS, flg);
			}
			if (!check)
				continue;

			// add filter after adjustment 20120912
			float originalTratio = snv.getTumor().getRatio();

			boolean initFilterFalse = originalTratio < KarkinosProp.mintumorratio;
			boolean highnormalAfterTCadjust = false;
			boolean lowtumorAfterTCadjust = false;
			boolean lowdepthafterTCadjust = false;

			boolean lowdepthafterTCadjustInfo = false;
			int numSupportRead = snv.getTumor().getAltCnt();
			int numSupportReadsNormal = snv.getNormal().getAltCnt();
			int normalTotal = snv.getNormal().getTotalcnt();
			boolean isindel = snv.getFlg() == PileUP.TumorINDEL;
			// Fisher test
			int[] normala = snv.getNormal().getRefAltCnt();
			int[] tumora = snv.getTumor().getRefAltCnt();
			// adjusted by tumor rate add 20120804

			double logn = snv.getNormal().getRefLogLikeHood();
			double logt = snv.getTumor().getMutateLogLikeHoodAmongMutation();

			//
			String chrom = snv.getChr();
			int pos = snv.getPos();

			if (pos == debugpos) {
				System.out.println("here");
			}

			FilterResult fr = new FilterResult();

			fr.setNormalVrcnt(numSupportReadsNormal);
			fr.setTumorVrcnt(numSupportRead);
			fr.setIndel(isindel);

			fr.setIndel(isindel);
			if (isindel) {
				boolean uncertainDelation = false;
				if (!snv.getTumor().isInsersion()) {
					//
					uncertainDelation = snv.getTumor().getDelmap().size() > 1;
				} else {
					uncertainDelation = snv.getTumor().getInsersionmap().size() > 1;
				}
				fr.setUncertainDelation(uncertainDelation);
			}

			// 1.mappability check
			float mappability = 1;
			if (checkMappability) {
				mappability = mappabilityanno.getMappability(chrom, pos);
			}
			// 2.entropy check
			String before10 = tgr.getGenomicSeq(chrom, pos - 10, pos - 1, true);
			String neighbor10 = tgr.getGenomicSeq(chrom, pos - 5, pos - 1, true)
					+ tgr.getGenomicSeq(chrom, pos + 1, pos + 5, true);

			// 20170317
			boolean indelinrepeat = false;
			int start = pos + 1;
			int end = pos + 10;
			String del = "";
			if (isindel) {
				try {
					if (!snv.getTumor().isInsersion()) {
						del = snv.getTumor().getIndelStr().split("\t")[0];
						start = pos + del.length() + 1;
						end = pos + del.length() + 10;

					}
				} catch (Exception ex) {

				}
			}

			String after10 = tgr.getGenomicSeq(chrom, start, end, true);
			if (del.length() >= 3) {
				if (after10.toUpperCase().contains(del)) {
					indelinrepeat = true;
				}
			}

			double s_b4 = Entropy.entropy(before10);
			double s_mid = Entropy.entropy(neighbor10);
			double s_after = Entropy.entropy(after10);
			int getminidx = getminidx(s_b4, s_mid, s_after);
			double seqEntropy = 0;
			String sseq = "";
			if (getminidx == 0) {
				seqEntropy = s_b4;
				sseq = before10;
			} else if ((getminidx == 1)) {
				seqEntropy = s_mid;
				sseq = neighbor10;
			} else if ((getminidx == 2)) {
				seqEntropy = s_after;
				sseq = after10;
			}

			// high isoform
			boolean highisoform = false;
			if (ge != null) {

				highisoform = ge.isFrequentIsoform(chrom, pos);
				fr.setHighisoform(highisoform);
			}

			char altTumor = snv.getTumor().getALT();
			int cntalt = countAltOrder(sseq, altTumor);
			if (cntalt >= 2) {
				fr.setIgnoreEntroptCheck(true);
			}

			char ref = snv.getNormal().getGenomeR();
			char alt = snv.getNormal().getALT();
			String key = ref + "to" + alt;
			double snpratio = nsnpC.getHetroSNPRatio(key);
			double adjustedLogn = logn;

			double[] afs = RatioUtils.getRatio(snv, ref, alt);

			fr.setTumorAF(afs[0]);
			fr.setNormalAF(afs[1]);
			fr.setNormaldepth(snv.getNormal().getTotalcnt());
			fr.setTumordepth(snv.getTumor().getTotalcnt());
			fr.setLogn(logn);
			fr.setLogt(logt);

			fr.setMapQ(snv.getTumor().getMQrms());
			fr.setDbSNPbean(snv.getDbSNPbean());
			fr.setMappability(mappability);
			fr.setSeqEntropy(seqEntropy);
			// bayesian scoreing to be impl
			// fr.setBresult(bfilter.scoring(snv));
			fr.setPhredbq(snv.getTumor().getPhred());
			snv.getTumor().getPhredQual();
			fr.setSecondAllele(snv.getTumor().getSecondALT());
			fr.setSecondAF(snv.getTumor().getSecondRatio());

			snv.setFilterResult(fr);

		}

	}

	public void filterAnnotation(DataSet dataset, ReadsSummary readsSummary, int ploidy) throws IOException {

		// BayesianFilter bfilter = new
		// BayesianFilter(dataset.getAnalyseDist());
		NormalSNPCounter nsnpC = new NormalSNPCounter(dataset);

		Map<String, Integer> snppos = new HashMap<String, Integer>();
		// float tumorContentsRatio = dataset.getAnalyseDist().getTumorratio();
		// fix 2013/07/24

		float tumorContentsRatio = dataset.getTumorRatio();

		int allcount = 0;
		int oxoGcnt = 0;
		int ffpecnt = 0;

		for (SNVHolder snv : dataset.getSnvlist()) {
			//
			//
			allcount++;
			boolean oxoGCand = SupportReadsCheck.oxoG(snv.getTumor().getGenomeR(), snv.getTumor().getALT());
			boolean ffpeCand = SupportReadsCheck.ffpe(snv.getTumor().getGenomeR(), snv.getTumor().getALT());
			if (oxoGCand) {
				oxoGcnt++;
			}
			if (ffpeCand) {
				ffpecnt++;
			}

		}
		float oxoGratio = (float) ((double) oxoGcnt / (double) allcount);
		float ffpeRatio = (float) ((double) ffpecnt / (double) allcount);

		System.out.println("oxoGratio =" + oxoGratio + " ffpe=" + ffpeRatio);

		for (SNVHolder snv : dataset.getSnvlist()) {

			 if(snv.getPos()==debugpos){
			 System.out.println("here");
			 }

			int flg = snv.getFlg();
			boolean check = (flg == PileUP.SomaticMutation || snv.getFlg() == PileUP.TumorINDEL);
			if (flg == PileUP.REGBOTH || flg == PileUP.NormalSNP || flg == PileUP.BothINDEL
					|| flg == PileUP.NormalINDEL) {
				String snpposS = snv.getChr() + "-" + snv.getPos();
				snppos.put(snpposS, flg);
			}
			if (!check)
				continue;

			// add filter after adjustment 20120912
			float originalTratio = snv.getTumor().getRatio();
			float adjustedTumorAllereFreq = CalcUtils.getTumorrateAdjustedRatio(snv, tumorContentsRatio);
			int adjustedtumortotalreads = Math.round(snv.getTumor().getTotalcnt() * tumorContentsRatio);

			boolean initFilterFalse = originalTratio < KarkinosProp.mintumorratio;
			boolean highnormalAfterTCadjust = false;
			boolean lowtumorAfterTCadjust = false;
			boolean lowdepthafterTCadjust = false;

			boolean lowdepthafterTCadjustInfo = false;

			// boolean lowdepth = snv.getTumor().getTotal() <
			// KarkinosProp.tcFilterdepth;
			// boolean lowfreq = adjustedTumorAllereFreq <
			// KarkinosProp.tcFilterfreq;
			int numSupportRead = snv.getTumor().getAltCnt();
			int numSupportReadsNormal = snv.getNormal().getAltCnt();
			int normalTotal = snv.getNormal().getTotalcnt();
			boolean isindel = snv.getFlg() == PileUP.TumorINDEL;

			String chrom = snv.getChr();
			int pos = snv.getPos();

			char altTumor = snv.getTumor().getALT();
			String neighbor60 = tgr.getGenomicSeq(chrom, pos - 30, pos + 30, true);
			boolean typicalSysErr = IlluminaSysError.checkTypicalError(snv.getTumor().getGenomeR(), altTumor,
					neighbor60, numSupportRead);

			if (initFilterFalse) {

				// read is resqued by low initial thres
				float nr = snv.getNormal().getRatio();
				float normalThresWithTumorContentsPenalty = KarkinosProp.maxnormalratio;
				if (tumorContentsRatio > 0 && tumorContentsRatio < 1) {
					normalThresWithTumorContentsPenalty = (KarkinosProp.maxnormalratio * tumorContentsRatio);
				}
				if (nr > normalThresWithTumorContentsPenalty) {
					highnormalAfterTCadjust = true;
				}
				if (numSupportRead < 10) {
					if (snv.getNormal().getAltCnt() > 0) {
						highnormalAfterTCadjust = true;
					}
				}

				int readspenalty = 1;
				// if (isindel) {
				// readspenalty = 1;
				// }
				if (numSupportRead < KarkinosProp.minsupportreads + readspenalty) {
					// lowdepthafterTCadjust = true;
					// mod 20130710 change to 2nd filter
					lowdepthafterTCadjustInfo = true;

				}
				if (numSupportRead == KarkinosProp.minsupportreads + readspenalty) {

					if (snv.getTumor().getRatio() < (KarkinosProp.min_initial_tumorratio + 0.02)) {
						lowdepthafterTCadjustInfo = true;
					}
				}

			}

			///////////////////////////////////
			//
			// Fisher test
			int[] normala = snv.getNormal().getRefAltCnt();
			int[] tumora = snv.getTumor().getRefAltCnt();
			// adjusted by tumor rate add 20120804
			int tumorref = Math.round(tumora[0] * tumorContentsRatio);
			double pvalFisher = Fisher.calcPValue(normala[0], normala[1], tumorref, tumora[1]);
			boolean fisherTestSignif = pvalFisher <= KarkinosProp.Fisher_Thres_For_SNV_Detection;
			snv.setFisherTestSignif(fisherTestSignif);
			snv.setPvalFisher(pvalFisher);

			double logn = snv.getNormal().getRefLogLikeHood();
			double logt = snv.getTumor().getMutateLogLikeHoodAmongMutation();

			//

			FilterResult fr = new FilterResult();

			fr.setNormalVrcnt(numSupportReadsNormal);
			fr.setTumorVrcnt(numSupportRead);
			fr.setIndel(isindel);

			fr.setIndel(isindel);
			if (isindel) {
				boolean uncertainDelation = false;
				if (!snv.getTumor().isInsersion()) {
					//
					uncertainDelation = snv.getTumor().getDelmap().size() > 1;
				} else {
					uncertainDelation = snv.getTumor().getInsersionmap().size() > 1;
				}
				fr.setUncertainDelation(uncertainDelation);
			}

			// 1.mappability check
			float mappability = 1;
			if (checkMappability) {
				mappability = mappabilityanno.getMappability(chrom, pos);
			}
			// 2.entropy check
			String before10 = tgr.getGenomicSeq(chrom, pos - 10, pos - 1, true);
			String neighbor10 = tgr.getGenomicSeq(chrom, pos - 5, pos - 1, true)
					+ tgr.getGenomicSeq(chrom, pos + 1, pos + 5, true);

			// 20170317 todai
			boolean indelinrepeat = false;
			int start = pos + 1;
			int end = pos + 10;
			
			int dellenn = 0;
			String del="";
			if (isindel) {
				try {
					if (!snv.getTumor().isInsersion()) {
						String dellens = snv.getTumor().getIndelStr().split("\t")[0];
						dellenn = Integer.parseInt(dellens);
						start = pos + dellenn;
						end = pos + dellenn + 10;
						del = tgr.getGenomicSeq(chrom, pos+1,pos+dellenn, true);
						
					}
				} catch (Exception ex) {

				}
			}

			String after10 = tgr.getGenomicSeq(chrom, start, end, true);
			if(pos==9710663){
				System.out.println("debug");
			}
			if (del.length() >= 3 && del.length() <=5) {
				if (after10.toUpperCase().contains(del.toUpperCase())) {
					indelinrepeat = true;
				}
			}

			double s_b4 = Entropy.entropy(before10);
			double s_mid = Entropy.entropy(neighbor10);
			double s_after = Entropy.entropy(after10);
			int getminidx = getminidx(s_b4, s_mid, s_after);
			double seqEntropy = 0;
			String sseq = "";
			if (getminidx == 0) {
				seqEntropy = s_b4;
				sseq = before10;
			} else if ((getminidx == 1)) {
				seqEntropy = s_mid;
				sseq = neighbor10;
			} else if ((getminidx == 2)) {
				seqEntropy = s_after;
				sseq = after10;
			}

			// high isoform
			boolean highisoform = false;
			if (ge != null) {

				highisoform = ge.isFrequentIsoform(chrom, pos);
				fr.setHighisoform(highisoform);
			}

			int cntalt = countAltOrder(sseq, altTumor);
			if (cntalt >= 2) {
				fr.setIgnoreEntroptCheck(true);
			}

			// 3.repeak check
			boolean repeat = containlowerCase(neighbor10);

			// 3.5 check typical illumina error

			// 4. support reads check
			SupportReadsCheckResult srr = srCheck.checkSupportReads(chrom, pos, snv.getTumor(), snv.getFlg(),
					dataset.getAnalyseDist().getTumorratio(), (snv.getCi().getVaridateVal() * 2), highisoform,
					normalTotal, snppos, adjustedTumorAllereFreq, snv.getNormal().getRatio(), oxoGratio, ffpeRatio);

			Set<Integer> supportReadsFlgs = srr.getFilter();
			fr.setPval4FiserDirectional(srr.getPval4directionCheck());
			if (repeat) {
				supportReadsFlgs.add(FilterResult.INFO_Maskedrepeat_BY_RepeatMasker);
			}
			// 20170317 todai
			if(indelinrepeat){
				
				if((normalTotal<snv.getTumor().getTotalcnt()*0.5) && (normalTotal<50) && snv.getTumor().getAltCnt()<=12){
				  supportReadsFlgs.add(FilterResult.Low_complexty);
				}
			}
			
			
			// not indel case, filter out
			// where normal low quality reads are rich
			if (!snv.getTumor().isIndel()) {
				if (snv.getNormal().getLowqualratio() > KarkinosProp.lowQualRatiothres) {
					//
					supportReadsFlgs.add(FilterResult.HighLowQualReads);
				}
			}

			char ref = snv.getNormal().getGenomeR();
			char alt = snv.getNormal().getALT();
			String key = ref + "to" + alt;
			// double adjustratio
			// =
			// nsnpC.getLogRefHetroSNPRatio(key,readsSummary.getNucCountRef(ref));

			// double snpratio = nsnpC.getHetroSNPRatio(key);
			// System.out.println("snpratio="+snpratio);
			// double adjustratio = Math.log10(snpratio / 0.1666);
			// double adjustratio = Math.log10(1-snpratio);
			// double adjustedLogn = adjustratio + logn;
			// 20130322 kill snp prior adjustment. no good effect
			double adjustedLogn = logn;

			// double ar = nsnpC.getHetroSNPRatioRemain(key);
			// double adjustedLogn = ar*logn;
			double[] afs = RatioUtils.getRatio(snv, ref, alt);
			fr.setSupportreadsBAlleleFeqquency(srr.getSupportreadsBAlleleFeqquency());
			fr.setRefreadsBAlleleFeqquency(srr.getRefreadsBAlleleFeqquency());
			fr.setTumorAF(afs[0]);
			fr.setNormalAF(afs[1]);
			fr.setInitFilterFalse(initFilterFalse);
			fr.setNumSupportRead(numSupportRead);
			fr.setNormaldepth(snv.getNormal().getTotalcnt());
			fr.setTumordepth(snv.getTumor().getTotalcnt());
			fr.setLogn(logn);
			fr.setLogt(logt);
			fr.setLogtAjusted(srr.getLogtAdjusted());
			fr.setMapQ(snv.getTumor().getMQrms());
			fr.setFisherTestSignif(fisherTestSignif);
			fr.setPvalFisher(pvalFisher);
			fr.setDbSNPbean(snv.getDbSNPbean());
			fr.setMappability(mappability);
			fr.setRepeat(repeat);
			fr.setSeqEntropy(seqEntropy);
			fr.setSupportReadsFlgs(supportReadsFlgs);

			// bayesian scoreing to be impl
			// fr.setBresult(bfilter.scoring(snv));
			fr.setPhredbq(snv.getTumor().getPhred());
			fr.setAdjustedTumorAllereFreq(adjustedTumorAllereFreq);
			fr.setAdjustedtumortotalreads(adjustedtumortotalreads);
			fr.setLowdepthafterTCadjustInfo(lowdepthafterTCadjustInfo);

			// TC adjust
			fr.setHighnormalAfterTCadjust(highnormalAfterTCadjust);
			// fr.setLowtumorAfterTCadjust(lowtumorAfterTCadjust);
			fr.setLowdepthafterTCadjust(lowdepthafterTCadjust);
			fr.setLognAjusted(adjustedLogn);

			snv.getTumor().getPhredQual();
			fr.setSecondAllele(snv.getTumor().getSecondALT());
			fr.setSecondAF(snv.getTumor().getSecondRatio());
			fr.setAfter10(after10);
			fr.setBefore10(before10);
			fr.setTypicalSysErr(typicalSysErr);

			snv.setFilterResult(fr);

		}

		// recaluculate psotrior probability of SNV
		SNVHighDepthCounter snvc = new SNVHighDepthCounter(dataset);
		na = new NoiseAnalysis();
		na.analysisNoiseRegion(dataset, ploidy);

		for (SNVHolder snv : dataset.getSnvlist()) {

			int flg = snv.getFlg();
			boolean check = (flg == PileUP.SomaticMutation || snv.getFlg() == PileUP.TumorINDEL);
			if (!check)
				continue;
			//
			char ref = snv.getNormal().getGenomeR();
			char alt = snv.getNormal().getALT();
			String key = ref + "to" + alt;
			key = GenotypeKeyUtils.aggrigateKeys(key);
			double snvratio = snvc.getSNVRatio(key);
			FilterResult fr = snv.getFilterResult();
			float originalTratio = snv.getTumor().getRatio();
			boolean initFilterFalse = originalTratio < KarkinosProp.mintumorratio;

			if (fr != null && !fr.isIndel()) {

				//
				// set depthAF matrix
				if (initFilterFalse) {
					boolean reject = false;
					if (na.reject(snv, tumorContentsRatio, fr.isIndel())) {
						reject = true;
						fr.setLowtumorAfterTCadjust(true);
					}
					System.out.println(
							snv.getChr() + "\t" + snv.getPos() + "\t" + snv.getTumor().getRatio() + "\t" + reject);

				}
				// 20130322 kill snv prior adjustment. no good effect
				// double adjustratio = Math.log10(snvratio);
				// //double adjustratio = Math.log10(1-snpratio);
				// double adjustedLogt = adjustratio + fr.getLogtAjusted();
				// fr.setLogtAjusted(adjustedLogt);
				fr.getPassFilterFlgForce();

			}

		}

	}

	public NoiseAnalysis getNa() {
		return na;
	}

	private int countAltOrder(String sseq, char altTumor) {

		int n = 0;
		try {
			int[] atgc = new int[4];
			String ref = "ATGC";
			for (char ch : sseq.toCharArray()) {

				ch = Character.toUpperCase(ch);
				int idx = ref.indexOf(ch);
				atgc[idx] = atgc[idx] + 1;
			}
			int idx = ref.indexOf(altTumor);
			int cntalt = atgc[idx];
			int m = 0;
			for (int cnt : atgc) {

				if (m == idx) {
					m++;
					continue;
				}
				if (cnt > cntalt) {
					n++;
				}
				m++;
			}
		} catch (Exception ex) {
		}
		return n;
	}

	private int getminidx(double d1, double d2, double d3) {

		if (d1 <= Math.min(d2, d3))
			return 0;
		if (d2 <= Math.min(d1, d3))
			return 1;
		if (d3 <= Math.min(d2, d1))
			return 2;
		return 0;
	}

	private boolean containlowerCase(String s) {

		for (char c : s.toCharArray()) {
			if (Character.isLowerCase(c)) {
				return true;
			}
		}
		return false;
	}

}
