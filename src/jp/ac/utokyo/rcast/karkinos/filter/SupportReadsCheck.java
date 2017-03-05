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

import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.CONTAIN_Reccurent_MISMATCH;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.INFO_AllelicInfoAvailable;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.INFO_SUPPORTED_BY_ONEDirection;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.LOW_PROPER;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.MutationAtSameCycle;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.NEARINDEL;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.NOISE_IN_NORMAL;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.READSENDSONLY;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.SUPPORTED_BY_ONEDirection;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.SUPPORTED_READSNG;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.TNQualityDiff;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.TooManyMismatchReads;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.noStrandSpecific;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.softClip;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.annotation.stats.Fisher;
import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPPool;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult.Counter;
import jp.ac.utokyo.rcast.karkinos.utils.OddRatioCalculator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class SupportReadsCheck extends ReadWriteBase {

	String bam;
	String normalbam;
	SAMFileReader bamr = null;
	SAMFileReader normalbamr = null;
	TwoBitGenomeReader tgr = null;
	DbSNPAnnotation dbAnno = null;

	public SupportReadsCheck(String normalbam, String bam,
			TwoBitGenomeReader tgr, DbSNPAnnotation dbAnno) {

		this.normalbam = normalbam;
		this.bam = bam;
		bamr = getReader(bam);
		normalbamr = getReader(normalbam);
		this.tgr = tgr;
		this.dbAnno = dbAnno;
	}

	public class SCounter implements java.io.Serializable {
		int n = 1;

		void inc() {
			n++;
		}
	}

	public static void main(String[] arg) {

		//
		// String bam =
		// "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.8/summary_cfDNA1/bam/Se-67-tumor25-Se-67_b-DNA_tumor_genome.bam";
		// String nbam =
		// "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/MT/MT30-2-MT30N/normal/MT30-2-MT30N_normal_genome.bam";
		// String bam =
		// "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/MT/MT30-2-MT30N/tumor/MT30-2-MT30N_tumor_genome.bam";
		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";

		String nbam = "/data3/users/yamamoto/exome/CRC/karkinos4.1.11/summary_CRC_all_samples/bam/CRC107_T-CRC107_N_normal_genome.bam";
		String bam = "/data3/users/yamamoto/exome/CRC/karkinos4.1.11/summary_CRC_all_samples/bam/CRC107_T-CRC107_N_tumor_genome.bam";
		String middelfile = "/GLUSTER_DIST/data/users/yamamoto/exome/CRC/karkinos2.0.3/CRC_107_T-CRC_107_N/sobj";

		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		DbSNPAnnotation dbAnno = null;

		SupportReadsCheck inst = new SupportReadsCheck(nbam, bam, tgr, dbAnno);
		//
		PileUPResult pileUPResult = new PileUPResult();
		pileUPResult.setGenomeRef('C');

		IndelInfo ii = new IndelInfo();
		// ii.indel = true;
		// ii.length = 17;
		// ii.cnt = 13;
		// // ii.insersion = "G";

		pileUPResult.setBaseAndQual('A', (byte) 50, 1, ii);
		Map<String, Integer> snppos = new HashMap<String, Integer>();
		try {
			inst.checkSupportReads("chr12", 25398284, pileUPResult,
					PileUP.SomaticMutation, 0.2f, PileUP.SomaticMutation,
					false, 10, snppos, 0.2f, 0.2f, 0f, 0f);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public SupportReadsCheckResult checkSupportReads(String chr, int pos,
			PileUPResult pileUPResult, int pileupFlg, float tumorratioss,
			double copynumber, boolean highisoform, int normalTotal,
			Map<String, Integer> snppos, float adjustedTumorAllereFreq,
			float normalAF, float oxoGratio, float ffpeRatio)
			throws IOException {

		// boolean debug = false;

		boolean debug = false;
		int debugpos = 7577114;
		if (pos == debugpos) {
			debug = true;
		}

		SupportReadsCheckResult ret = new SupportReadsCheckResult();
		Set<Integer> filter = new LinkedHashSet<Integer>();
		ret.setFilter(filter);

		CloseableIterator<SAMRecord> ite = bamr.query(chr, pos, pos, false);

		List<SAMRecord> supportreads = new ArrayList<SAMRecord>();
		List<SAMRecord> referencereads = new ArrayList<SAMRecord>();

		List<SAMRecord> allreads = new ArrayList<SAMRecord>();

		int rightelse = 0;
		int leftelse = 0;
		int rightsupport = 0;
		int leftsupport = 0;
		boolean pairedend = false;
		int proper = 0;
		boolean isindel = pileUPResult.isIndel();

		List<BaseandQData> pileuplist = new ArrayList<BaseandQData>();
		Map<Integer, SCounter> cycleCheck = new HashMap<Integer, SCounter>();

		int f1r2support = 0;
		int f2r1support = 0;

		int f1r2else = 0;
		int f2r1else = 0;

		while (ite.hasNext()) {

			SAMRecord sam = ite.next();
			// System.out.println(sam.getReadString());
			allreads.add(sam);
			int[] reta = containTargetMutation(pos, pileUPResult, sam,
					pileupFlg);
			int mutationidx = reta[1];
			boolean samContatinMutation = (reta[0] == 0);
			// check log likehood t
			int readslen = sam.getReadLength();
			if ((mutationidx >= 0) && (mutationidx < readslen)) {
				char ch = (char) sam.getReadBases()[mutationidx];
				byte qual = sam.getBaseQualities()[mutationidx];
				if ((int) qual < KarkinosProp.minPhredQualForEach) {
					continue;
				}

				pileuplist.add(new BaseandQData(ch, qual));
				//
			}
			int flg = sam.getFlags();
			// System.out.println(flg);

			if (samContatinMutation) {

				// System.out.println(sam.getReadName());
				// if (sam.getFirstOfPairFlag()) {
				if (sam.getReadNegativeStrandFlag()) {
					leftsupport++;

					if (sam.getFirstOfPairFlag()) {
						f2r1support++;
					} else {
						f1r2support++;
					}

				} else {
					rightsupport++;

					if (sam.getFirstOfPairFlag()) {
						f1r2support++;
					} else {
						f2r1support++;
					}

				}

				if (sam.getReadPairedFlag()) {
					pairedend = true;
					if (sam.getProperPairFlag()) {
						proper++;
					}
				}

				// } else {
				// if (sam.getReadNegativeStrandFlag()) {
				// rightsupport++;
				// } else {
				// leftsupport++;
				// }
				// }

				SCounter cyclecounter = null;
				if (cycleCheck.containsKey(mutationidx)) {
					cyclecounter = cycleCheck.get(mutationidx);
					cyclecounter.inc();
				} else {
					cyclecounter = new SCounter();
					cycleCheck.put(mutationidx, cyclecounter);
				}

				supportreads.add(sam);

			} else {

				// if (sam.getFirstOfPairFlag()) {
				if (sam.getReadNegativeStrandFlag()) {
					leftelse++;

					if (sam.getFirstOfPairFlag()) {
						f2r1else++;
					} else {
						f1r2else++;
					}

				} else {
					rightelse++;

					if (sam.getFirstOfPairFlag()) {
						f1r2else++;
					} else {
						f2r1else++;
					}

				}
				// } else {
				// if (sam.getReadNegativeStrandFlag()) {
				// rightelse++;
				// } else {
				// leftelse++;
				// }
				// }

			}

		}
		ite.close();

		// add near indel check around SNV 2012.08.03
		if (!isindel) {
			try {
				boolean neighborIndelExist = neighborIndelExist(allreads, pos);
				if (neighborIndelExist) {
					filter.add(NEARINDEL);
				}
				if (!isindel) {
					double r = (double) pileUPResult.indelcnt()
							/ (double) supportreads.size();
					//
					if (r > 0.1) {
						filter.add(NEARINDEL);
					}
				}

			} catch (Exception ex) {
			}
		}
		float adjustedLogt = 100;
		float ntQualityDiff = 100;
		if (!isindel) {
			// 20130319 use fixed copynunber 1
			float fixedcn = 1.0f;
			adjustedLogt = getAdjuatedLogt(pileUPResult.getGenomeR(),
					pileUPResult.getALT(), pileuplist, tumorratioss, fixedcn,
					normalAF);
			ret.setLogtAdjusted(adjustedLogt);

			ntQualityDiff = getNTQDiff(pileUPResult.getGenomeR(),
					pileUPResult.getALT(), pileuplist);

			ret.setNtQualityDiff(ntQualityDiff);

			float qdiffthres = KarkinosProp.TNQdiff;
			boolean highsupport = supportreads.size() > 10;
			boolean midsupport = supportreads.size() > 6;
			if (highsupport) {
				qdiffthres = qdiffthres + 4;
			} else if (midsupport) {
				qdiffthres = qdiffthres + 2;
			}

			if (ntQualityDiff > qdiffthres) {
				filter.add(TNQualityDiff);
			}
		}
		//
		//

		String mtypeparent = pileUPResult.getGenomeR() + "-"
				+ pileUPResult.getALT();

		if (supportreads.size() == 0) {
			filter.add(SUPPORTED_READSNG);
			return ret;
		}

		double fisherP2 = 1;
		if (f1r2support == 0 || f2r1support == 0) {
			// return SUPPORTED_BY_ONEDirection;
			filter.add(INFO_SUPPORTED_BY_ONEDirection);
			fisherP2 = Fisher.calcPValue(f1r2else, f2r1else, f1r2support,
					f2r1support);

			boolean f1r2 = (f2r1support == 0);
			boolean oxoGCand = oxoG(pileUPResult.getGenomeR(),
					pileUPResult.getALT(), f1r2);
			boolean ffpeCand = ffpe(pileUPResult.getGenomeR(),
					pileUPResult.getALT(), f1r2);

			double thres = KarkinosProp.Fisher_Thres_For_Reads_Direction;

			if (oxoGCand) {

				if (f1r2support <= 6 && f2r1support <= 6) {

					if (oxoGratio > 0.1) {
						thres = KarkinosProp.Fisher_Thres_For_Reads_Direction2;
					}

				}

				// if (oxoGratio > 0.2) {
				// thres = KarkinosProp.Fisher_Thres_For_Reads_Direction3;
				// }

				// System.out.println("oxoG" + pileUPResult.getGenomeR() +" "
				// +pileUPResult.getALT() + "P=" + fisherP2+" " +f1r2else+" "
				// +f2r1else+" "+ f1r2support+" "+f2r1support);

			}
			if (ffpeCand) {

				if (f1r2support <= 8 && f2r1support <= 8) {
					
					if (ffpeRatio > 0.1) {
						thres = KarkinosProp.Fisher_Thres_For_Reads_Direction2;
					}
					if (ffpeRatio > 0.2) {
						thres = KarkinosProp.Fisher_Thres_For_Reads_Direction3;
					}
				}
				// if (f1r2support >= 6 || f2r1support >= 6) {
				// thres = 1;
				// }

				// System.out.println("ffpe" + pileUPResult.getGenomeR() +" "
				// +pileUPResult.getALT() + "P=" + fisherP2+" " +f1r2else+" "
				// +f2r1else+" "+ f1r2support+" "+f2r1support);

			}
			//

			boolean fisherSignif = (fisherP2 <= thres);
			if (fisherSignif) {
				filter.add(SUPPORTED_BY_ONEDirection);
			}
		}

		double fisherP = 1;

		if (leftsupport == 0 || rightsupport == 0) {
			// return SUPPORTED_BY_ONEDirection;
			filter.add(INFO_SUPPORTED_BY_ONEDirection);
			fisherP = Fisher.calcPValue(rightelse, leftelse, rightsupport,
					leftsupport);
			boolean fisherSignif = fisherP <= KarkinosProp.Fisher_Thres_For_Reads_Direction;
			if (fisherSignif) {
				filter.add(SUPPORTED_BY_ONEDirection);
			}

			// if (cycleCheck.size() <= 1) {
			// if (!pileUPResult.isIndel()) {
			// filter.add(MutationAtSameCycle);
			// }
			// }
			// mod 2013/07/18
			float ssRatio = sameCycleRatio(cycleCheck);
			if (ssRatio >= 0.75) {
				if (!pileUPResult.isIndel()) {
					filter.add(MutationAtSameCycle);
				}
			}

		}
		ret.setPval4directionCheck((float) fisherP);

		float mismatchRate = getMismatchRate(supportreads, pos);
		if (mismatchRate > KarkinosProp.minMisMatchRate) {
			filter.add(TooManyMismatchReads);
		}

		float softCliprate = getSoftClipRate(supportreads);
		if (softCliprate > 0.8) {
			filter.add(softClip);
		}

		for (SAMRecord sam : allreads) {

			if (!supportreads.contains(sam)) {
				referencereads.add(sam);
			}
		}

		int size = supportreads.size();
		int rangestart = supportreads.get(0).getAlignmentStart();
		int rangeend = supportreads.get(size - 1).getAlignmentEnd();
		//
		List[] la = getPairList(bamr, supportreads, referencereads, chr,
				rangestart, rangeend, KarkinosProp.pairReadsMaxdist);
		List<SAMRecord> supportreadsPair = la[0];
		List<SAMRecord> otherPair = la[1];

		int sizep = supportreadsPair.size();
		if (sizep == 0) {
			filter.add(SUPPORTED_READSNG);
			return ret;
		}
		if (pairedend) {

			double properratio = (double) proper / (double) sizep;
			if (properratio < 0.3) {
				filter.add(LOW_PROPER);
			}

		}
		// normal bam check
		if (isindel) {

			boolean normalUnClear = normalCheck(chr, pos, pileUPResult);
			if (normalUnClear) {
				filter.add(NOISE_IN_NORMAL);
			}

		}

		int rangestartp = supportreadsPair.get(0).getAlignmentStart();
		int rangeendp = supportreadsPair.get(sizep - 1).getAlignmentEnd();

		// DbSNPBean dbSNPBeanPos = dbAnno.lookup(chr, pos);

		int recnt = 0;
		int diffcountneighbor = 0;
		int diffrecnt = 0;
		//
		//
		for (int n = rangestartp; n <= rangeendp; n++) {

			// debugpos = 26671523;
			// //System.out.println(n);
			// if (n == debugpos) {
			// System.out.println(pos);
			// }

			// System.out.println(n);
			PileUPResult mutationsupport = PileUPPool.borrowObject();
			PileUPResult refsupport = PileUPPool.borrowObject();
			char genomerefn = tgr.getGenomeNuc(chr, n, true);
			genomerefn = Character.toUpperCase(genomerefn);
			mutationsupport.setGenomeRef(genomerefn);
			refsupport.setGenomeRef(genomerefn);
			// /

			int readsendmutation = setPileup(chr, n, genomerefn,
					mutationsupport, supportreadsPair, isindel);
			if (n == pos) {
				if (readsendmutation >= 0 && readsendmutation <= 3) {
					filter.add(READSENDSONLY);
				}
				// do not pileup at target pos
				continue;
			}
			setPileup(chr, n, genomerefn, refsupport, otherPair, isindel);

			boolean mismatchng = mutationsupport.isDiff()
					&& !refsupport.isDiff();
			float tumorratio = getTratio(mutationsupport, refsupport);
			int tumorcnt = getTcnt(mutationsupport);
			boolean enoughmismatch = ((tumorratio) >= (pileUPResult.getRatio() * 0.5))
					&& tumorcnt > 2;
			DbSNPBean dbSNPBean = null;
			if (dbAnno != null) {
				dbSNPBean = dbAnno.lookup(chr, n);
			}
			boolean dbSNP = isSNP(dbSNPBean);
			if (mutationsupport.isDiff()) {
				if (enoughmismatch && !dbSNP) {
					diffcountneighbor++;
				}
			}

			if (mismatchng) {

				// System.out.println(n+"\t"+mutationsupport.getRefAltCnt()[0]+"\t"+mutationsupport.getRefAltCnt()[1]);

				boolean enoughrefreads = refsupport.getTotalcnt() > 2;
				// DbSNPBean dbSNPBean = dbAnno.lookup(chr, n);

				String poss = chr + "-" + n;
				boolean snpPosE = snppos.containsKey(poss);

				if (enoughmismatch && enoughrefreads) {

					if (!dbSNP && !snpPosE) {

						recnt++;
						// only count different substitution type
						String mtype = genomerefn + "-"
								+ mutationsupport.getALT();
						boolean sametype = mtype.equalsIgnoreCase(mtypeparent)
								&& (Math.abs(pos - n) <= 5);
						if (!sametype) {
							diffrecnt++;
						}
					}
				}
				// check neighbor SNP allele frequency
				//
				boolean bothsnp = false;
				if (snpPosE) {
					bothsnp = (snppos.get(poss) == PileUP.REGBOTH || snppos
							.get(poss) == PileUP.BothINDEL);
				}
				if (dbSNPBean != null && bothsnp) {

					filter.add(INFO_AllelicInfoAvailable);
					boolean candidatehetro = (adjustedTumorAllereFreq <= 0.5);
					float supportreadsBAlleleFeqquency = mutationsupport
							.getRatio();
					float refreadsBAlleleFeqquency = refsupport.getRatio();
					ret.setRefreadsBAlleleFeqquency(refreadsBAlleleFeqquency);
					ret.setSupportreadsBAlleleFeqquency(supportreadsBAlleleFeqquency);
					//
					if (candidatehetro && hetro(mutationsupport)
							&& hetro(refsupport)) {
						filter.add(noStrandSpecific);
					}

				}

			}

			//

			//
			PileUPPool.returnObject(mutationsupport);
			PileUPPool.returnObject(refsupport);

		}

		// boolean regrec = false;
		// if (recnt >= 2) {
		// regrec = true;
		// } else if (recnt == 1) {
		// // boolean dbSNP = (dbSNPBeanPos != null &&
		// !dbSNPBeanPos.isCosmic());
		// // boolean lowdepthN = (normalTotal <=
		// KarkinosProp.low_normal_depth_thres);
		// // if(dbSNP || lowdepthN || highisoform){
		// // regrec = true;
		// // }
		// }
		if (recnt >= 2) {
			filter.add(CONTAIN_Reccurent_MISMATCH);
		}
		if (diffrecnt >= 1) {
			filter.add(CONTAIN_Reccurent_MISMATCH);
		}
		if (diffcountneighbor >= 8) {
			filter.add(CONTAIN_Reccurent_MISMATCH);
		}
		// change value 2013.07.03
		// if (recnt >= 3) {
		// filter.add(CONTAIN_Reccurent_MISMATCH);
		// }
		// if (diffrecnt >= 2) {
		// filter.add(CONTAIN_Reccurent_MISMATCH);
		// }

		return ret;

	}

	public static boolean ffpe(char genomeR, char alt, boolean f1r2) {

		// if (genomeR == 'C' && alt == 'T' && f1r2) {
		// return true;
		// }
		// if (genomeR == 'G' && alt == 'A' && !f1r2) {
		// return true;
		// }
		if (genomeR == 'C' && alt == 'T') {
			return true;
		}
		if (genomeR == 'G' && alt == 'A') {
			return true;
		}
		return false;
	}

	public static boolean oxoG(char genomeR, char alt, boolean f1r2) {

		if (genomeR == 'C' && alt == 'A' && f1r2) {
			return true;
		}
		if (genomeR == 'G' && alt == 'T' && !f1r2) {
			return true;
		}
		return false;

	}

	public static boolean oxoG(char genomeR, char alt) {
		// TODO Auto-generated method stub
		return oxoG(genomeR, alt, true) || oxoG(genomeR, alt, false);
	}

	public static boolean ffpe(char genomeR, char alt) {

		return ffpe(genomeR, alt, true) || ffpe(genomeR, alt, false);
	}

	//
	private boolean normalCheck(String chr, int pos, PileUPResult pileUPResult)
			throws IOException {

		boolean ret = false;

		int s = 0;
		int e = 1;
		int len = 0;
		String insstr = "";
		if (pileUPResult.isInsersion()) {

			Iterator<String> ite = pileUPResult.getInsersionmap().keySet()
					.iterator();
			while (ite.hasNext()) {
				String ss = ite.next();
				if (ss.length() > len) {
					len = ss.length();
					insstr = ss;
				}
			}
			s = -len;
			e = len;
		} else {

			Iterator<Integer> ite = pileUPResult.getDelmap().keySet()
					.iterator();

			while (ite.hasNext()) {
				int lenm = ite.next();
				if (lenm > len) {
					len = lenm;
				}
			}
			s = 0;
			e = len;

		}

		CloseableIterator<SAMRecord> ite = normalbamr.query(chr, pos + s, pos
				+ e, false);

		List<SAMRecord> list = new ArrayList<SAMRecord>();
		List<SAMRecord> mismatchL = new ArrayList<SAMRecord>();
		//
		while (ite.hasNext()) {

			SAMRecord sam = ite.next();
			list.add(sam);
			Integer nm = sam.getIntegerAttribute("NM");
			if (nm != null && nm >= 3) {
				mismatchL.add(sam);
			}

		}
		ite.close();

		// try mismatch resque
		int resquedbase = 0;
		for (SAMRecord sam : mismatchL) {

			int mismatchb4 = sam.getIntegerAttribute("NM");
			int mismatchafter = mismatchb4;

			try {
				mismatchafter = mismatchAfterIndel(chr, pos, sam, pileUPResult,
						mismatchb4, len, insstr);
			} catch (Exception ex) {

			}

			int diff = mismatchb4 - mismatchafter;
			if (diff > 0) {
				resquedbase = resquedbase + diff;
			}

		}
		if (resquedbase >= 8) {
			return true;
		}
		//
		int cnt = 0;
		for (int n = pos + s; n <= pos + e; n++) {

			//
			PileUPResult pu = PileUPPool.borrowObject();
			// private int setPileup(String chr, int pos, char genomeR,
			// PileUPResult ret,
			// List<SAMRecord> reads, boolean isindel) {

			char genomerefn = tgr.getGenomeNuc(chr, n, true);
			pu.setGenomeRef(genomerefn);
			setPileup(chr, n, genomerefn, pu, list, false);
			if (pu.getRatio() > 0.05) {

				boolean dbSNP = false;
				if (dbAnno != null) {
					DbSNPBean dbSNPBean = dbAnno.lookup(chr, n);
					dbSNP = isSNP(dbSNPBean);
				}
				if (!dbSNP) {
					cnt++;
				}
			}

		}
		// false if more than 5% freq
		return cnt >= 2;

	}

	private int mismatchAfterIndel(String chr, int pos, SAMRecord sam,
			PileUPResult pileUPResult, int mismatchb4, int len, String insstr)
			throws IOException {

		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.DELETION)) {
				return mismatchb4;
			}

			if (ce.getOperator().equals(CigarOperator.INSERTION)) {
				return mismatchb4;
			}
		}

		int mism = 0;
		if (pileUPResult.isInsersion()) {

			int start = sam.getAlignmentStart();
			int end = sam.getAlignmentEnd();

			for (int n = start; n <= end; n++) {

				if (n < pos) {

					int cidx = getCharIdx(n, sam, new IndelInfo());
					char genomerefn = tgr.getGenomeNuc(chr, n, true);
					if ((cidx >= 0) && (cidx < sam.getReadLength())) {
						char c = sam.getReadString().charAt(cidx);
						if (c != genomerefn) {

							mism++;
						}
					}

				} else if (n > pos) {

					int cidx = getCharIdx(n + len, sam, new IndelInfo());
					char genomerefn = tgr.getGenomeNuc(chr, n, true);
					if ((cidx >= 0) && (cidx < sam.getReadLength())) {
						char c = sam.getReadString().charAt(cidx);
						if (c != genomerefn) {

							mism++;
						}
					}

				} else {

					// n
					int cidx = getCharIdx(n, sam, new IndelInfo());
					int id = 0;
					for (int m = cidx; m < cidx + len; m++) {
						if ((m >= 0) && (m < sam.getReadLength())) {
							char c = sam.getReadString().charAt(m);
							char c2 = insstr.charAt(id);
							if (c != c2) {
								mism++;
							}
						}
						id++;
					}
				}

			}

		} else {

			if (pos <= sam.getAlignmentStart()
					&& sam.getAlignmentStart() <= pos + len) {

				int start = sam.getAlignmentStart() - len;
				int end = sam.getAlignmentEnd();
				// System.out.println(sam.getReadName());
				// System.out.println(sam.getReadString());

				for (int n = start; n <= end; n++) {

					if (n > pos + len) {

						char genomerefn = tgr.getGenomeNuc(chr, n, true);
						int cidx = getCharIdx(n, sam, new IndelInfo());
						if ((cidx >= 0) && (cidx < sam.getReadLength())) {
							char c = sam.getReadString().charAt(cidx);
							if (c != genomerefn) {

								mism++;
							}
						}

					} else {

						if (n < pos + len) {
							continue;
						}
						char genomerefn = tgr.getGenomeNuc(chr, n, true);
						int cidx = getCharIdx(n + len, sam, new IndelInfo());
						if ((cidx >= 0) && (cidx < sam.getReadLength())) {
							char c = sam.getReadString().charAt(cidx);
							if (c != genomerefn) {
								mism++;
							}
						}

					}
				}

			} else {

				int start = sam.getAlignmentStart();
				int end = sam.getAlignmentEnd() + len;
				// System.out.println(sam.getReadName());
				// System.out.println(sam.getReadString());

				for (int n = start; n <= end; n++) {

					if (n < pos) {

						char genomerefn = tgr.getGenomeNuc(chr, n, true);
						int cidx = getCharIdx(n, sam, new IndelInfo());
						if ((cidx >= 0) && (cidx < sam.getReadLength())) {
							char c = sam.getReadString().charAt(cidx);
							if (c != genomerefn) {

								mism++;
							}
						}
					} else {

						if (n < pos + len) {
							continue;
						}
						char genomerefn = tgr.getGenomeNuc(chr, n, true);
						int cidx = getCharIdx(n - len, sam, new IndelInfo());
						if ((cidx >= 0) && (cidx < sam.getReadLength())) {
							char c = sam.getReadString().charAt(cidx);
							if (c != genomerefn) {
								mism++;
							}
						}

					}

				}

			}

		}
		System.out.println(sam.getAlignmentStart() + " "
				+ sam.getAlignmentEnd() + " " + "mismatchb4=" + mismatchb4
				+ " mismatch" + mism + " " + sam.getReadName());
		return mism;

	}

	private boolean isSNP(DbSNPBean dbSNPBean) {

		if (dbSNPBean == null) {
			return false;
		} else {
			if (dbSNPBean.isValid()) {
				return true;
			}
		}
		return false;
	}

	private float getNTQDiff(char genomeR, char alt,
			List<BaseandQData> pileuplist) {

		SummaryStatistics ssref = new SummaryStatistics();
		SummaryStatistics ssalt = new SummaryStatistics();
		for (BaseandQData data : pileuplist) {

			//
			if (data.base == genomeR) {
				ssref.addValue(data.qual);
			} else if (data.base == alt) {
				ssalt.addValue(data.qual);
			}

		}

		double refave = 10;
		double altave = 10;
		if (ssref.getN() >= 10) {
			refave = ssref.getMean();
		}
		altave = ssalt.getMean();
		double diff = refave - altave;
		return (float) diff;
	}

	private float sameCycleRatio(Map<Integer, SCounter> cycleCheck) {

		Iterator<Integer> ite = cycleCheck.keySet().iterator();
		int max = 0;
		int total = 0;
		while (ite.hasNext()) {

			//
			SCounter sc = cycleCheck.get(ite.next());
			total = total + sc.n;

			if (sc.n > max) {
				max = sc.n;
			}

		}
		double r = (double) max / (double) total;
		return (float) r;
	}

	private boolean hetro(PileUPResult ppr) {
		float r = ppr.getRatio();
		return (r > 0.4 && r < 0.6) && (ppr.getTotal() >= 15);
	}

	private List<SAMRecord>[] getPairList(SAMFileReader bamr2,
			List<SAMRecord> supportreads, List<SAMRecord> referencereads,
			String chr, int start, int end, int margin) {

		CloseableIterator<SAMRecord> ite = null;
		List[] rary = new List[2];
		rary[0] = supportreads;
		rary[1] = referencereads;
		try {
			Set<String> nameset = new HashSet<String>();
			for (SAMRecord sr : supportreads) {
				nameset.add(sr.getReadName());
			}

			List<SAMRecord> supportreadsPair = new ArrayList<SAMRecord>();
			List<SAMRecord> otherPair = new ArrayList<SAMRecord>();
			int s = start - margin;
			if (s < 1)
				s = 1;
			int e = end + margin;
			ite = bamr.query(chr, s, e, false);
			while (ite.hasNext()) {

				SAMRecord sam = ite.next();
				if (nameset.contains(sam.getReadName())) {
					//
					supportreadsPair.add(sam);
				} else {
					//
					otherPair.add(sam);
				}

			}
			rary[0] = supportreadsPair;
			rary[1] = otherPair;

			if (supportreadsPair.size() <= supportreads.size()) {
				// should not happen
				rary[0] = supportreads;
			}
			if (otherPair.size() <= referencereads.size()) {
				// should not happen
				rary[1] = referencereads;
			}

		} catch (Exception ex) {

		} finally {
			if (ite != null) {
				try {
					ite.close();
				} catch (Exception ex) {

				}
			}
		}

		return rary;

	}

	private String key(SAMRecord sam) {
		boolean pflg = false;
		if (sam.getReadPairedFlag()) {
			pflg = sam.getFirstOfPairFlag();
		}
		return sam.getReadName() + pflg;
	}

	private boolean neighborIndelExist(List<SAMRecord> allreads, int pos) {

		int mergin = KarkinosProp.nearindelbt;
		int total = 0;
		int indelc = 0;
		for (SAMRecord sam : allreads) {

			int[] indelIndex = getIndelIdx(sam);
			total++;
			if (indelIndex == null) {
				continue;
			} else {
				int indelposs = sam.getAlignmentStart() + indelIndex[0];
				int indelpose = sam.getAlignmentStart() + indelIndex[0]
						+ indelIndex[1];

				if (Math.abs(indelposs - pos) == 0
						|| (Math.abs(indelpose - pos) == 0)) {
					continue;
				}

				boolean b1 = (Math.abs(indelposs - pos) < mergin);
				boolean b2 = (Math.abs(indelpose - pos) < mergin);
				if (b1 || b2) {
					indelc++;
				}
			}

		}
		double ratio = (double) indelc / (double) total;
		return ratio > 0.05;// more than 5%
	}

	// culculate adjusted tumor ratio
	private float getAdjuatedLogt(char genomeR, char alt,
			List<BaseandQData> pileuplist, double tumorratio,
			double copynumber, double normalAF) {
		//

		return OddRatioCalculator.getAdjuatedLogt(genomeR, alt, pileuplist,
				tumorratio, copynumber, normalAF);
	}

	private float getSoftClipRate(List<SAMRecord> supportreads) {

		try {
			int scc = 0;
			for (SAMRecord sam : supportreads) {

				boolean sc = false;
				for (CigarElement ce : sam.getCigar().getCigarElements()) {
					if (ce.getOperator().equals(CigarOperator.SOFT_CLIP)) {
						sc = true;
						break;
					}
				}
				if (sc) {
					scc++;
				}

			}
			double d = (double) scc / (double) supportreads.size();
			return (float) d;
		} catch (Exception e) {

		}
		return 0f;
	}

	private float getMismatchRate(List<SAMRecord> supportreads, int pos)
			throws IOException {

		int totalbase = 0;
		int mismatchbase = 0;
		IndelInfo indelinfo = new IndelInfo();
		//
		for (SAMRecord sam : supportreads) {

			int start = sam.getAlignmentStart();
			int end = sam.getAlignmentEnd();
			String reads = sam.getReadString();
			if (end == 0) {
				end = start + sam.getReadLength();
			}
			for (int n = start; n <= end; n++) {

				if (n == pos)
					continue;
				int cidx = getCharIdx(n, sam, indelinfo);
				if ((cidx >= 0) && (cidx < sam.getReadLength())) {
					byte qual = sam.getBaseQualities()[cidx];
					if ((int) qual < KarkinosProp.minPhredQualForEach) {
						continue;
					}
					char read = reads.charAt(cidx);
					char ref = tgr
							.getGenomeNuc(sam.getReferenceName(), n, true);

					totalbase++;
					if (notEqualChar(read, ref)) {

						mismatchbase++;
					}
				}
			}

		}
		return (float) ((double) mismatchbase / (double) totalbase);

	}

	private boolean notEqualChar(char c1, char c2) {

		c1 = Character.toUpperCase(c1);
		c2 = Character.toUpperCase(c2);
		return c1 != c2;
	}

	private int getTcnt(PileUPResult mutationsupport) {
		char alt = mutationsupport.getALT();
		int idx = PileUPResult.seqALL.indexOf(alt);
		int[] sq = mutationsupport.getSeqCounter();
		int diffcnt = sq[idx];
		return diffcnt;
	}

	private float getTratio(PileUPResult mutationsupport,
			PileUPResult refsupport) {

		int total = mutationsupport.getTotalcnt() + refsupport.getTotalcnt();
		int diffcnt = getTcnt(mutationsupport);

		return (float) ((double) diffcnt / (double) total);
	}

	private int setPileup(String chr, int pos, char genomeR, PileUPResult ret,
			List<SAMRecord> reads, boolean isindel) {

		int maxreadsend = -1;
		boolean diff = false;
		IndelInfo indelinfo = new IndelInfo();
		indelinfo.indel = isindel;
		indelinfo.refpos = -1;
		int depth = 0;

		for (SAMRecord sam : reads) {

			indelinfo.clear();
			indelinfo.refpos = -1;
			int seqidx = getCharIdx(pos, sam, indelinfo);
			char ch = 0;
			byte qual = 0;
			int readslen = sam.getReadLength();
			int readsend = 0;

			if ((seqidx >= 0) && (seqidx < readslen)) {
				ch = (char) sam.getReadBases()[seqidx];
				qual = sam.getBaseQualities()[seqidx];
				readsend = Math.min(readslen - seqidx, seqidx);
				if (isindel) {
					readsend = Math.min(readslen - indelinfo.refpos,
							indelinfo.refpos);
				}
				int mq = sam.getMappingQuality();
				ret.setBaseAndQual(ch, qual, mq, indelinfo);
				if (maxreadsend < readsend) {
					maxreadsend = readsend;
				}
				if (ch != genomeR) {
					diff = true;
				}
				depth++;
			}

		}
		ret.setDiff(diff);

		return maxreadsend;
	}

	// return index
	private int[] containTargetMutation(int pos, PileUPResult pileUPResult,
			SAMRecord sam, int pileupFlg) {

		boolean ret = false;
		IndelInfo indelInfo = new IndelInfo();
		int seqidx = getCharIdx(pos, sam, indelInfo);
		// tumor somatic mutation
		char ch = 0;
		byte qual = 0;
		if ((seqidx >= 0) && (seqidx < sam.getReadLength())) {
			ch = (char) sam.getReadBases()[seqidx];
			qual = sam.getBaseQualities()[seqidx];
		}

		if (pileupFlg == PileUP.TumorINDEL) {

			if (indelInfo.isIndel()) {
				String insersion = indelInfo.getInsersion();
				if (insersion == null) {
					// delation
					Map<Integer, Counter> m = pileUPResult.getDelmap();
					if (m == null) {
						ret = false;
					} else {
						ret = m.containsKey(indelInfo.length);
					}
				} else {
					// insersion
					Map<String, Counter> m = pileUPResult.getInsersionmap();
					if (m == null) {
						ret = false;
					} else {
						ret = m.containsKey(insersion);
					}
				}

			}

		} else {

			boolean targetEqual = (ch == pileUPResult.getALT());
			ret = targetEqual;

		}
		int[] retary = new int[3];
		if (seqidx < 0)
			seqidx = 0;
		// if(ret){
		// if(seqidx<0)seqidx=0;
		// return seqidx;
		// }else{
		// return -1;
		// }
		retary[0] = ret ? 0 : 1;
		retary[1] = seqidx;
		return retary;
	}

	private static int getCharIdx(int pos, SAMRecord sam, IndelInfo indelinfo) {

		indelinfo.indel = false;
		int start = sam.getAlignmentStart();
		int relpos = pos - start;
		if (relpos == 0)
			return 0;

		int readidx = 0;
		int refidx = 0;

		List<CigarElement> list = sam.getCigar().getCigarElements();
		int l = 0;
		for (CigarElement ce : list) {

			int len = ce.getLength();
			if (len == sam.getReadLength()) {
				return relpos;
			}

			if (ce.getOperator().consumesReferenceBases()
					&& ce.getOperator().consumesReadBases()) {

				if (relpos <= refidx + len) {

					int readidxr = readidx + (relpos - refidx);
					// check if insersion exsist in next cigar
					if (relpos == refidx + len) {
						if (l + 1 < list.size()) {
							CigarElement nextcigar = list.get(l + 1);
							if (nextcigar.getOperator() == CigarOperator.INSERTION) {
								indelinfo.indel = true;
								indelinfo.length = nextcigar.getLength();
								indelinfo.insersion = substring(
										sam.getReadString(), readidxr, readidxr
												+ indelinfo.length);
								indelinfo.refpos = refidx + len;
							} else if (nextcigar.getOperator() == CigarOperator.DELETION) {
								indelinfo.indel = true;
								indelinfo.length = nextcigar.getLength();
								indelinfo.refpos = refidx + len;
							}
						}
					}
					return readidxr;
				}
				refidx += len;
				readidx += len;

			} else if (ce.getOperator().consumesReferenceBases()) {

				if (ce.getOperator() == CigarOperator.DELETION) {
					if (relpos == refidx + len) {
						indelinfo.indel = true;
						indelinfo.length = len;
						indelinfo.refpos = refidx + len;
						return -1;
					} else if (relpos < refidx + len) {
						return -1;
					}
				}
				if (ce.getOperator() == CigarOperator.N) {

					if (relpos < refidx + len) {
						return -1;
					}
				}
				refidx += len;
			} else {

				readidx += len;

			}
			l++;
		}
		return readidx;

	}

	private int[] getIndelIdx(SAMRecord sam) {
		List<CigarElement> list = sam.getCigar().getCigarElements();
		int relpos = 0;
		for (CigarElement ce : list) {

			int len = ce.getLength();
			if (len == sam.getReadLength()) {
				return null;
			}

			if ((ce.getOperator() == CigarOperator.DELETION)
					|| (ce.getOperator() == CigarOperator.INSERTION)) {

				int dellen = ce.getOperator() == CigarOperator.DELETION ? len
						: 0;
				return new int[] { relpos, dellen };

			} else if (ce.getOperator().consumesReadBases()) {

				relpos += len;

			}
		}

		return null;

	}

	private static String substring(String str, int s, int e) {

		if (e >= str.length()) {
			return str.substring(s);
		} else {
			return str.substring(s, e);
		}

	}

}
