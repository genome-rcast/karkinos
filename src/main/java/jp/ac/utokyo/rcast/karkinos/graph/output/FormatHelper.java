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
package jp.ac.utokyo.rcast.karkinos.graph.output;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.filter.FilterResult;
import jp.ac.utokyo.rcast.karkinos.filter.FilterSNP;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class FormatHelper {

	// bw.append("##fileformat=VCFv4.1"+ "\n");
	//
	// bw.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total
	// Depth\">");
	// bw.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele
	// Frequency adjusted by tumor ratio\">");
	// bw.append("##INFO=<ID=AFO,Number=A,Type=Float,Description=\"Allele
	// Frequency original\">");
	// bw.append("##INFO=<ID=CN,Number=A,Type=Float,Description=\"Copy
	// Number\">");
	// bw.append("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP
	// membership\">");
	// "##INFO=<ID=pval,Number=0,Type=Float,Description=\"Fisher test pval\">"
	// bw.append("##FILTER=<ID=qf,Description=\"Quality below threshold\">");
	// bw.append("##FILTER=<ID=snp,Description=\"dbSNP snp\">");
	// bw.append("##FILTER=<ID=ma,Description=\"Low mappability < 0.5\">");
	// bw.append("##FILTER=<ID=entropy,Description=\"Low complexty\">");
	// bw.append("##FILTER=<ID=srd,Description=\"Support reads have only one
	// direction\">");
	// bw.append("##FILTER=<ID=srm,Description=\"Only support reads have another
	// mismatch\">");
	// bw.append("##FILTER=<ID=sre,Description=\"mutation at reads ends\">");
	//
	// bw.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	// bw.append("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype
	// Quality\">");
	// bw.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read
	// Depth\">");
	//
	// bw.append("#CHROM POS ID REF ALT QUAL FILTER INFO FILTER2");
	//

	public static String getVCFLine(SNVHolder snv, TwoBitGenomeReader tgr, float tumorRratio, NoiseAnalysis na)
			throws IOException {

		int flg = snv.getFlg();
		// boolean reg = (flg == PileUP.SomaticMutation || flg ==
		// PileUP.TumorINDEL);
		// if (!reg) {
		// return null;
		// }

		if (flg == PileUP.SomaticMutation) {
			String[] data = FormatHelper.getVCFCol(snv, tumorRratio, na);
			return FormatHelper.tabDelimated(data);
		} else if (flg == PileUP.TumorINDEL) {
			String[] data = FormatHelper.getVCFColIndel(snv, tgr, tumorRratio, na);
			return FormatHelper.tabDelimated(data);
		}
		return null;

	}

	public static String[] getLine(SNVHolder snv, TwoBitGenomeReader tgr, GeneExons ge) {

		String[] data = new String[16];
		// chr
		data[0] = snv.getChr().replace("chr", "");
		// pos
		data[1] = String.valueOf(snv.getPos());
		String freq = "";
		if (snv.getDbSNPbean() != null) {
			// id
			data[2] = snv.getDbSNPbean().getInfo();
			if (snv.getDbSNPbean().getMode() == DbSNPAnnotation.MODE1000g) {
				//
				freq = String.valueOf(snv.getDbSNPbean().getFreq());
			}
		} else {
			data[2] = ".";
		}
		PileUPResult pir = snv.getNormal();
		char genome = pir.getGenomeR();
		// ref
		data[3] = String.valueOf(Character.toUpperCase(genome));
		// alt
		data[4] = String.valueOf(Character.toUpperCase(pir.getALT()));
		// qual
		data[5] = String.valueOf((pir.getPhred() / 10));
		// refcnt;
		data[6] = String.valueOf(pir.getRefCnt());
		// altcnt;
		data[7] = String.valueOf(pir.getAltCnt());
		// freq;
		data[8] = String.valueOf(pir.getRatio());

		PileUPResult pirt = snv.getTumor();
		// alt
		data[9] = String.valueOf(Character.toUpperCase(pir.getALT()));
		// qual
		data[10] = String.valueOf((pirt.getPhred() / 10));
		// refcnt;
		data[11] = String.valueOf(pirt.getRefCnt());
		// altcnt;
		data[12] = String.valueOf(pirt.getAltCnt());
		// freq;
		data[13] = String.valueOf(pirt.getRatio());

		// info
		data[14] = getInfoStr4Normal(snv, ge, freq);

		// tn ratio
		data[15] = String.valueOf(getInfoTN_AFratio(snv));
		//

		return data;

	}

	public static String getVCFLineAllDiff(SNVHolder snv, TwoBitGenomeReader tgr, GeneExons ge) {

		int flg = snv.getFlg();
		try {

			if (flg == PileUP.REGBOTH || flg == PileUP.NormalSNP || flg == PileUP.NONSignif
					|| flg == PileUP.SomaticMutation) {

				String[] data = FormatHelper.getVCFCol4Normal(snv, ge);
				if (validate(data)) {
					return FormatHelper.tabDelimated(data);
				}

			} else if (flg == PileUP.BothINDEL || flg == PileUP.NormalINDEL || flg == PileUP.TumorINDEL) {
				String[] data = FormatHelper.getVCFColIndel4Normal(snv, tgr, ge);
				if (validate(data)) {
					return FormatHelper.tabDelimated(data);
				}
			}

		} catch (Exception ex) {

		}
		return null;
	}

	private static boolean validate(String[] data) {

		return !data[3].equals(data[4]);

	}

	public static String getVCFLine4SNP(SNVHolder snv, TwoBitGenomeReader tgr, GeneExons ge) throws IOException {

		int flg = snv.getFlg();

		if (flg == PileUP.REGBOTH || flg == PileUP.NormalSNP) {
			String[] data = FormatHelper.getVCFCol4Normal(snv, ge);
			return FormatHelper.tabDelimated(data);
		} else if (flg == PileUP.BothINDEL || flg == PileUP.NormalINDEL) {
			String[] data = FormatHelper.getVCFColIndel4Normal(snv, tgr, ge);
			return FormatHelper.tabDelimated(data);
		}
		return null;

	}

	public static String[] getVCFCol4Normal(SNVHolder snv, GeneExons ge) {

		String[] data = new String[9];
		// chr
		data[0] = snv.getChr().replace("chr", "");
		// pos
		data[1] = String.valueOf(snv.getPos());
		String freq = "";
		if (snv.getDbSNPbean() != null) {
			// id
			data[2] = snv.getDbSNPbean().getInfo();
			if (snv.getDbSNPbean().getMode() == DbSNPAnnotation.MODE1000g) {
				//
				freq = String.valueOf(snv.getDbSNPbean().getFreq());
			}
		} else {
			data[2] = ".";
		}
		PileUPResult pir = snv.getNormal();
		char genome = pir.getGenomeR();
		// ref
		data[3] = String.valueOf(Character.toUpperCase(genome));
		// alt
		data[4] = String.valueOf(Character.toUpperCase(pir.getALT()));
		// qual
		data[5] = String.valueOf((pir.getPhred() / 10));

		data[6] = FilterSNP.filter(pir);
		// info
		data[7] = getInfoStr4Normal(snv, ge, freq);
		// tn ratio
		data[8] = String.valueOf(getInfoTN_AFratio(snv));
		return data;
	}

	public static String[] getVCFCol(SNVHolder snv, float tumorRratio, NoiseAnalysis na) {

		String[] data = new String[13];
		// chr
		data[0] = snv.getChr().replace("chr", "");
		// pos
		data[1] = String.valueOf(snv.getPos());
		if (snv.getDbSNPbean() != null) {
			// id
			data[2] = snv.getDbSNPbean().getInfo();
		} else {
			data[2] = ".";
		}
		PileUPResult pir = snv.getTumor();
		char genome = pir.getGenomeR();
		// ref
		data[3] = String.valueOf(Character.toUpperCase(genome));
		// alt
		data[4] = String.valueOf(Character.toUpperCase(pir.getALT()));
		// qual
		data[5] = String.valueOf(pir.getPhred());
		// filter
		FilterResult fr = snv.getFilterResult();
		Set<Integer> flgs = null;
		String fs = "PASS";
		String b4 = "";
		String after = "";
		if (fr != null) {

			if (!fr.isPassFilter()) {
				flgs = fr.getPassFilterFlg();
				String fsa = fr.getFilterStr(flgs);
				if (fsa.length() > 1) {
					fs = fsa;
				} else {
					fs = "PASS";
				}

			}
			b4 = fr.getBefore10();
			after = fr.getAfter10();
		}
		data[6] = fs;
		// info
		String f2 = getFilter2Str(snv, fs);
		float score = getScore(snv, fs);
		data[7] = getInfoStr(snv, tumorRratio, na, score);
		data[8] = f2;
		// normal
		data[9] = getAlleleStr(snv.getNormal());
		// tumor
		data[10] = getAlleleStr(snv.getTumor());
		//
		data[11] = b4;
		data[12] = after;

		return data;
	}

	private static float getScore(SNVHolder snv, String fs) {

		float score = 0;
		int errorcount = snv.getFilterResult().getPassFilterFlg().size();
		if (errorcount > 2) {
			score = 0;
		}
		if (errorcount == 2) {
			score = 0.2f;
			score = gradient(0.002, score, snv.getTumor().getAltCnt());
		}
		if (errorcount == 1) {
			score = 0.5f;
			score = gradient(0.002, score, snv.getTumor().getAltCnt());
		}

		if (fs.equals("PASS")) {
			score = 0.8f;
			score = gradient(0.01, score, snv.getTumor().getAltCnt());
		}
		String ps = getFilter2Str(snv, fs);
		if (ps.equals("PASS")) {
			score = 0.95f;
			score = gradient(0.01, score, snv.getTumor().getAltCnt());
		}
		return score;
	}

	private static float gradient(double r, float score, int i) {

		score = (float) (score + (r * i));
		if (score > 1)
			score = 1;
		return score;
	}

	private static String getAlleleStr(PileUPResult pl) {

		int[] refalt = pl.getRefAltCnt();
		return refalt[0] + ";" + refalt[1];
	}

	public static String[] getVCFColIndel(SNVHolder snv, TwoBitGenomeReader tgr, float tumorRratio, NoiseAnalysis na)
			throws IOException {
		String[] data = new String[13];
		data[0] = snv.getChr().replace("chr", "");
		// data[1] = String.valueOf(snv.getPos());
		if (snv.getDbSNPbean() != null) {
			data[2] = snv.getDbSNPbean().getInfo();
		} else {
			data[2] = ".";
		}
		PileUPResult pir = snv.getTumor();
		StringBuffer genome = new StringBuffer();
		StringBuffer alt = new StringBuffer();
		int adjustedpos = snv.getPos();

		if (pir.isInsersion()) {
			//
			genome.append(tgr.getGenomicSeq(snv.getChr(), snv.getPos() - 1, snv.getPos() - 1, true));
			String[] insersions = pir.getIndelStr().split("\t");
			for (String s : insersions) {
				if (alt.length() > 0) {
					alt.append(",");
				}
				alt.append(genome);
				alt.append(s);
			}
			data[1] = String.valueOf(adjustedpos);
			//
		} else {
			//

			String[] delations = pir.getIndelStr().split("\t");
			for (String s : delations) {
				if (genome.length() > 0) {
					genome.append(",");
				}
				int n = 0;
				try {
					n = Integer.parseInt(s);
				} catch (Exception ex) {

				}
				if (alt.length() == 0) {
					alt.append(tgr.getGenomicSeq(snv.getChr(), snv.getPos() - 1, snv.getPos() - 1, true));
				}
				String ss = tgr.getGenomicSeq(snv.getChr(), snv.getPos() - 1, snv.getPos() - 1 + n, true);
				genome.append(ss);

				data[1] = String.valueOf(adjustedpos);
			}

		}
		// ref
		data[3] = String.valueOf(genome.toString().toUpperCase());
		// alt
		data[4] = String.valueOf(alt.toString().toUpperCase());
		// qual
		data[5] = String.valueOf(".");
		// filter
		FilterResult fr = snv.getFilterResult();

		Set<Integer> flgs = null;
		String fs = "PASS";
		String b4 = "";
		String after = "";
		if (fr != null) {

			if (!fr.isPassFilter()) {

				flgs = fr.getPassFilterFlg();
				String fsa = fr.getFilterStr(flgs);
				if (fsa.length() > 1) {
					fs = fsa;
				}
			} else {
				fs = "PASS";
			}
			b4 = fr.getBefore10();
			after = fr.getAfter10();
		}
		data[6] = fs;
		// info
		float score = getScore(snv, fs);
		data[7] = getInfoStr(snv, tumorRratio, na, score);
		data[8] = getFilter2Str(snv, fs);
		// normal
		data[9] = getAlleleStr(snv.getNormal());
		// tumor
		data[10] = getAlleleStr(snv.getTumor());
		//
		data[11] = b4;
		data[12] = after;
		return data;
	}

	public static String[] getVCFColIndel4Normal(SNVHolder snv, TwoBitGenomeReader tgr, GeneExons ge)
			throws IOException {
		String[] data = new String[9];
		data[0] = snv.getChr().replace("chr", "");
		data[1] = String.valueOf(snv.getPos());
		String freq = "";
		if (snv.getDbSNPbean() != null) {
			// id
			data[2] = snv.getDbSNPbean().getInfo();
			if (snv.getDbSNPbean().getMode() == DbSNPAnnotation.MODE1000g) {
				//
				freq = String.valueOf(snv.getDbSNPbean().getFreq());
			}
		} else {
			data[2] = ".";
		}
		PileUPResult pir = snv.getNormal();
		StringBuffer genome = new StringBuffer();
		StringBuffer alt = new StringBuffer();

		if (pir.isInsersion()) {
			//
			genome.append(tgr.getGenomicSeq(snv.getChr(), snv.getPos() - 1, snv.getPos() - 1, true));
			String[] insersions = pir.getIndelStr().split("\t");
			for (String s : insersions) {
				if (alt.length() > 0) {
					alt.append(",");
				}
				alt.append(genome);
				alt.append(s);
			}

			//
		} else {
			//

			String[] delations = pir.getIndelStr().split("\t");
			for (String s : delations) {
				if (genome.length() > 0) {
					genome.append(",");
				}
				int n = Integer.parseInt(s);
				if (alt.length() == 0) {
					alt.append(tgr.getGenomicSeq(snv.getChr(), snv.getPos() - 1, snv.getPos() - 1, true));
				}
				String ss = tgr.getGenomicSeq(snv.getChr(), snv.getPos() - 1, snv.getPos() + Math.abs(n - 1), true);
				genome.append(ss);

			}

		}
		// ref
		data[3] = String.valueOf(genome.toString().toUpperCase());
		// alt
		data[4] = String.valueOf(alt.toString().toUpperCase());
		// qual
		data[5] = String.valueOf(".");
		// filter
		FilterResult fr = snv.getFilterResult();

		Set<Integer> flgs = null;
		String fs = FilterSNP.filterIndel(pir);
		data[6] = fs;
		// info
		data[7] = getInfoStr4Normal(snv, ge, freq);
		// data[8] = fs;
		// tn ratio
		data[8] = String.valueOf(getInfoTN_AFratio(snv));
		return data;
	}

	private static String getFilter2Str(SNVHolder snv, String fs) {

		boolean reff = snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_LOW_refOddsRatio);
		boolean trf = snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_LOW_tumorOddsRatio);
		boolean low = snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_adjustAlleleFreq);

		boolean lowd = snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_adjustLowdepth);

		boolean minsupport = snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_minimumSupportReads);

		boolean ffpe = snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_ffpe);

		boolean oxoG = snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_oxoG);

		if ((snv.getDbSNPbean() != null) && snv.getDbSNPbean().isCosmicvalid()) {
			// if filter 1 pass and cosmic then pass
			if (fs.equals("PASS")) {
				return fs;
			}
		}

		//
		if (reff || trf || low || lowd || minsupport || ffpe || oxoG) {
			StringBuffer sb = new StringBuffer();
			if (!fs.equals("PASS")) {
				sb.append(fs);
			}
			if (reff) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append("highref");
			}
			if (trf) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append("lowtumor");
			}
			if (low) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append("lowallelefreq");
			}
			if (lowd) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append("lowadsdepth");
			}
			if (minsupport) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append("lowsupportreads");
			}
			if (ffpe) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append("ffpe");
			}
			if (oxoG) {
				if (sb.length() > 0) {
					sb.append(",");
				}
				sb.append("oxoG");
			}
			return sb.toString();
		} else {
			return fs;
		}

	}

	// bw.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total
	// Depth\">");
	// bw.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele
	// Frequency adjusted by tumor ratio\">");
	// bw.append("##INFO=<ID=AFO,Number=A,Type=Float,Description=\"Allele
	// Frequency original\">");
	// bw.append("##INFO=<ID=CN,Number=A,Type=Float,Description=\"Copy
	// Number\">");
	// bw.append("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP
	// membership\">");
	private static String getInfoStr(SNVHolder snv, float tumorRratio, NoiseAnalysis na, float score) {

		StringBuffer is = new StringBuffer();
		//
		PileUPResult pir = snv.getTumor();

		is.append("DP=" + pir.getTotalcnt());
		float adjustedratio = CalcUtils.getTumorrateAdjustedRatio(snv, tumorRratio);
		is.append(",pvSNV=" + na.getPval(adjustedratio, (snv.getTumor().getTotalcnt() * tumorRratio)));
		is.append(",AF=" + adjustedratio);
		is.append(",AFO=" + pir.getRatio());
		is.append(",CN=" + (float) (snv.getCi().getVaridateVal() * 2));
		is.append(",ACN=" + (float) (snv.getCi().getAafreq()));
		is.append(",BCN=" + (float) (snv.getCi().getBafreq()));
		is.append(",BQ=" + (float) pir.getBQrms());
		is.append(",MQ=" + (float) pir.getMQrms());
		is.append(",S=" + (float) snv.getFilterResult().getSeqEntropy());
		is.append(",MP=" + (float) snv.getFilterResult().getMappability());
		is.append(",pD=" + (float) snv.getFilterResult().getPval4FiserDirectional());
		is.append(",ND=" + (float) snv.getNormal().getTotalcnt());
		is.append(",NR=" + (float) snv.getNormal().getRatio());
		is.append(",logn=" + (float) snv.getFilterResult().getLogn());
		is.append(",logt=" + (float) snv.getFilterResult().getLogt());
		is.append(",lognAdjuated=" + (float) snv.getFilterResult().getLognAjusted());
		is.append(",logtAdjuated=" + (float) snv.getFilterResult().getLogtAjusted());
		is.append(",pval=" + (float) snv.getPvalFisher());
		is.append(",SA=" + snv.getTumor().getSecondALT());
		is.append(",SAF=" + (float) snv.getTumor().getSecondRatio());
		is.append(",SCORE=" + score);

		if (snv.getDbSNPbean() != null) {
			if (snv.getDbSNPbean().isCosmic()) {
				// is.append(";cosmic");
				is.append(",cosmc=" + snv.getDbSNPbean().getCosmiccount());
			} else {
				is.append(",DBv=" + snv.getDbSNPbean().getVaridationStr());
			}
		}
		if (snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_AllelicInfoAvailable)) {

			is.append(",mBAF=" + (float) snv.getFilterResult().getSupportreadsBAlleleFeqquency());
			is.append(",sBAF=" + (float) snv.getFilterResult().getRefreadsBAlleleFeqquency());
			is.append(";AIA");

		}

		if (snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_SUPPORTED_BY_ONEDirection)) {
			is.append(";OD");
		}
		if (snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_Maskedrepeat_BY_RepeatMasker)) {
			is.append(";RP");
		}
		if (snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_LOW_refOddsRatio)) {
			is.append(";LowRefFilter");
		}
		if (snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_LOW_tumorOddsRatio)) {
			is.append(";LowTumorFilter");
		}
		if (snv.getDbSNPbean() != null) {
			if (snv.getDbSNPbean().isCosmic()) {
				is.append(";cosmic");
				if (snv.getDbSNPbean().isCosmicvalid()) {
					is.append(";cosm_v");
				}
			} else {
				is.append(";DB");
			}
		}
		if (snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_INIT_Filter_FAIL)) {
			is.append(";initfail");
		}
		if (snv.getFilterResult().getInfoFlg().contains(FilterResult.INFO_SNP_flg)) {
			is.append(";validateSNP");
		}
		return is.toString();

	}

	private static float getInfoTN_AFratio(SNVHolder snv) {

		float r = 1f;

		try {
			float rn = snv.getNormal().getRatio();
			float rt = snv.getTumor().getRatio();
			r = (float) ((double) rt / (double) rn);
		} catch (Exception ex) {
		}

		return r;
	}

	private static String getInfoStr4Normal(SNVHolder snv, GeneExons ge, String freq) {

		StringBuffer is = new StringBuffer();
		//
		PileUPResult pir = snv.getNormal();
		is.append("DP=" + pir.getTotalcnt());
		is.append(",AF=" + pir.getRatio());
		is.append(",BQ=" + (float) pir.getBQrms());
		is.append(",MQ=" + (float) pir.getMQrms());
		float tr = 0f;
		try {
			tr = snv.getTumor().getRatio();
		} catch (Exception ex) {

		}
		is.append(",AFT=" + tr);

		String refseq = ge.getGeneId(snv.getChr(), snv.getPos());
		if (refseq != null) {
			is.append(",NM=" + refseq);
			if (ge.onCDS(snv.getChr(), snv.getPos())) {
				is.append(",onCDS=" + true);
			}
			if (!freq.equals("")) {
				is.append(",1000freq=" + freq);
			}
		}
		try {
			is.append(",ACN=" + (float) (snv.getCi().getAafreq()));
			is.append(",BCN=" + (float) (snv.getCi().getBafreq()));
		} catch (Exception ex) {
		}
		return is.toString();

	}

	public static String tabDelimated(String[] sa) {

		if (sa == null)
			return "";
		StringBuffer sb = new StringBuffer();
		for (String s : sa) {
			if (s == null)
				break;
			if (sb.length() > 0) {
				sb.append("\t");
				sb.append(s);
			} else {
				sb.append(s);
			}
		}
		return sb.toString();
	}

	public static String getAnnoverInputLine(SNVHolder snv, TwoBitGenomeReader tgr, float tumorRratio) {
		int flg = snv.getFlg();

		if (flg == PileUP.SomaticMutation) {
			String[] data = FormatHelper.getAnnoverInputCol(snv, tumorRratio);
			return FormatHelper.tabDelimated(data);
		} else if (flg == PileUP.TumorINDEL) {
			String[] data = FormatHelper.getAnnoverInputColIndel(snv, tgr, tumorRratio);
			return FormatHelper.tabDelimated(data);
		}
		return null;
	}

	private static String[] getAnnoverInputCol(SNVHolder snv, float tumorRratio) {
		String[] data = new String[9];
		// chr
		data[0] = snv.getChr().replace("chr", "");
		// pos
		data[1] = String.valueOf(snv.getPos());
		data[2] = String.valueOf(snv.getPos());
		PileUPResult pir = snv.getTumor();
		char genome = pir.getGenomeR();
		// ref
		data[3] = String.valueOf(Character.toUpperCase(genome));
		// alt
		data[4] = String.valueOf(Character.toUpperCase(pir.getALT()));
		// qual
		data[5] = "unknown";
		// qual
		data[6] = String.valueOf(pir.getPhred());
		data[7] = String.valueOf(pir.getTotalcnt());
		data[8] = String.valueOf((float) pir.getMQrms());
		return data;
	}

	private static String[] getAnnoverInputColIndel(SNVHolder snv, TwoBitGenomeReader tgr, float tumorRratio) {
		String[] data = new String[9];
		// chr
		data[0] = snv.getChr().replace("chr", "");
		// pos

		PileUPResult pir = snv.getTumor();
		StringBuffer genome = new StringBuffer();
		if (pir.isInsersion()) {
			//
			data[1] = String.valueOf(snv.getPos());
			data[2] = String.valueOf(snv.getPos());
			data[3] = "-";
			String[] insersions = pir.getIndelStr().split("\t");
			String ins = insersions[0].toUpperCase();
			ins = ins.replaceAll("N", "A");
			data[4] = String.valueOf(ins);
			//
		} else {
			//

			String[] delations = pir.getIndelStr().split("\t");
			int n = 0;
			for (String s : delations) {
				if (genome.length() > 0) {
					genome.append(",");
				}
				n = Integer.parseInt(s);
				break;
			}
			String ss = "";
			try {
				ss = tgr.getGenomicSeq(snv.getChr(), snv.getPos(), snv.getPos() + Math.abs(n - 1), true);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			data[1] = String.valueOf(snv.getPos());
			data[2] = String.valueOf(snv.getPos() + (n - 1));
			data[3] = ss.toUpperCase();
			data[4] = "-";

		}

		// qual
		data[5] = "unknown";
		// qual
		data[6] = String.valueOf(".");
		// depth
		data[7] = String.valueOf(pir.getTotalcnt());
		data[8] = String.valueOf((float) pir.getMQrms());
		return data;
	}

}
