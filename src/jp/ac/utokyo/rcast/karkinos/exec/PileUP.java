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

import static jp.ac.utokyo.rcast.karkinos.exec.PileUPResult.seqALL;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;

import jp.ac.utokyo.rcast.karkinos.filter.DefinedSites;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.utils.Interval;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class PileUP implements java.io.Serializable {
	static boolean debug = false;
	public static int debugpos = 179604901;

	public static void pileup(Interval iv, DataSet dataset,
			List<SamHolder> normalList, List<SamHolder> tumorList,
			TwoBitGenomeReader tgr, ReadsSummary readsSummary, DefinedSites ds)
			throws IOException {

		//
		// LinkedList<SAMRecord> normalSubList = new LinkedList<SAMRecord>();
		// LinkedList<SAMRecord> tumorSubList = new LinkedList<SAMRecord>();

		SortedMap<Integer, CapInterval> sm = dataset.getCh()
				.getIntersectCapinterval(iv);
		if (sm == null)
			return;
		Iterator<Integer> ite = sm.keySet().iterator();
		// int normalSublistidx = 0;
		// int tumorSublistidx = 0;
		int mergin = KarkinosProp.baitmergin;

		int lastpileup = 0;
		while (ite.hasNext()) {

			//
			int key = ite.next();
			CapInterval ci = sm.get(key);
			int start = ci.start - mergin;
			if (lastpileup + 1 > start) {
				start = lastpileup + 1;
			}
			int end = ci.end + mergin;
			String chr = ci.chr;
			int[] indexes = new int[2];
			
			boolean definedpos = false;
			
			
			for (int n = start; n <= end; n++) {

				if(ds!=null){
					definedpos = ds.contains(n);
				}
				
				boolean ontag = ((n >= ci.start) && (n <= ci.end));
				
				if (iv.getStart() <= n && n <= iv.getEnd()) {

					indexes = _pileup(chr, n, dataset, normalList, tumorList,
							tgr, readsSummary, indexes, false,definedpos,ontag);

//					if (n == debugpos) {
//						System.out.println("here");
//						indexes = _pileup(chr, n, dataset, normalList,
//								tumorList, tgr, readsSummary, indexes, true);
//					}

					lastpileup = n;
				}

			}

		}

	}

	// private static int setSublist(int pos, LinkedList<SAMRecord> subList,
	// int idx, List<SAMRecord> list) {
	//
	// // remove
	// SAMRecord sam = null;
	// while (!subList.isEmpty() && !include(pos, (sam = subList.getFirst()))) {
	// subList.removeFirst();
	// }
	//
	// // add and return idx
	// int n = idx;
	// if (n < 0)
	// n = 0;
	// boolean onceAdd = false;
	// for (; n < list.size(); n++) {
	// SAMRecord samt = list.get(n);
	// if (include(pos, samt)) {
	// subList.addLast(samt);
	// onceAdd = true;
	// } else if (!onceAdd) {
	// continue;
	// } else {
	// n--;
	// break;
	// }
	// }
	// return n;
	// }

	// private static TreeMap<Integer, Integer> getListIdx(List<SAMRecord>
	// samList) {
	//
	// TreeMap<Integer, Integer> map = new TreeMap<Integer, Integer>();
	// int cnt = 0;
	// for(SAMRecord sam:samList){
	// int start = sam.getAlignmentStart();
	// if(cnt%1000==0){
	// map.put(start, cnt);
	// }
	// cnt++;
	// }
	// return map;
	// }

	private static int samplingFirstReads(List<SAMRecord> tumorList) {

		if (tumorList != null && !tumorList.isEmpty()) {
			return tumorList.get(0).getReadLength();
		}
		return 100;
	}

	// ////
	public static final int NONSignif = 0;
	public static final int NormalSNP = 1;
	public static final int SomaticMutation = 2;
	public static final int REGBOTH = 3;
	public static final int NormalINDEL = 4;
	public static final int TumorINDEL = 5;
	public static final int BothINDEL = 6;
	
	public static final int DISIGNATE = 7;

	private static String flgtoStr(int flg) {

		if (flg == 1) {
			return "NormalSNP";
		} else if (flg == 2) {
			return "Somatic";
		} else if (flg == 3) {
			return "REGBOTH";
		} else if (flg == 4) {
			return "NormalIndel";
		} else if (flg == 5) {
			return "SomaticIndel";
		} else if (flg == 6) {
			return "IndelOTH";
		}

		return "";
	}

	private static int[] _pileup(String chr, int pos, DataSet dataset,
			List<SamHolder> normalList, List<SamHolder> tumorList,
			TwoBitGenomeReader tgr, ReadsSummary readsSummary, int[] indexes,
			boolean debug, boolean definedpos, boolean ontag) throws IOException {

		int[] ret = new int[2];

		char genomeR = tgr.getGenomeNuc(chr, pos, true);

		boolean repeatmask = Character.isLowerCase(genomeR);
		genomeR = Character.toUpperCase(genomeR);

		PileUPResult normal = PileUPPool.borrowObject();
		int[] nret = setPileUP(normal, chr, pos, genomeR, normalList,
				indexes[0]);
		int ndepth = nret[0];
		int nidx = nret[1];
		//int nontarget = nret[2];
		int nontarget = ontag?1:0;
		
		ret[0] = nidx;
		readsSummary.setNormalDepth(chr, pos, ndepth, nontarget);

		PileUPResult tumor = PileUPPool.borrowObject();
		int[] tret = setPileUP(tumor, chr, pos, genomeR, tumorList, indexes[1]);
		int tdepth = tret[0];
		int tidx = tret[1];
		//int tontarget = nret[2];
		
		ret[1] = tidx;
		readsSummary.setTumorDepth(chr, pos, tdepth, nontarget);
		readsSummary.setNucCount(genomeR);
		
		boolean targetSeq = dataset.targetSeq;
		boolean indelreg = false;
		if (tumor.indel || normal.indel) {

			PileUPResult tumor2 = tumor.getIndelCopy();
			tumor2.setDiff(true);
			PileUPResult normal2 = normal.getIndelCopy();
			normal2.setDiff(true);
			int flg2 = checkReg(normal2, tumor2,targetSeq);
			
			if(definedpos&&flg2 == NONSignif){
				flg2 = DISIGNATE;
			}			
			
			if (flg2 != NONSignif) {
				SNVHolder snvHolder = new SNVHolder();
				snvHolder.setChr(chr);
				snvHolder.setPos(pos);
				snvHolder.setNormal(normal2);
				snvHolder.setTumor(tumor2);
				int iflg = NormalINDEL;
				if (tumor.indel && normal.indel) {
					iflg = BothINDEL;
				} else if (tumor.indel) {
					iflg = TumorINDEL;
				}

				snvHolder.setFlg(iflg);
				dataset.addSNVHolder(snvHolder);
				indelreg = true;
				// System.out.println(tumor2.getRatio());

				System.out.println("indel" + "\t" + chr + "\t" + pos + "\t"
						+ tumor2.getIndelStr() + "\t" + normal2.getRatio()
						+ "\t" + tumor2.getIndelStr() + "\t"
						+ tumor2.getRatio() + "\t");

			}
		}

		if (indelreg == false) {
			
			int flg = checkReg(normal, tumor,targetSeq);
			if(definedpos&&flg == NONSignif){
				flg = DISIGNATE;
			}		
			
			if (flg == NONSignif) {



				PileUPPool.returnObject(normal);
				PileUPPool.returnObject(tumor);
				readsSummary.setNucCountRef(genomeR);

			} else {

				SNVHolder snvHolder = new SNVHolder();
				snvHolder.setChr(chr);
				snvHolder.setPos(pos);
				snvHolder.setNormal(normal);
				snvHolder.setTumor(tumor);
				snvHolder.setFlg(flg);
				dataset.addSNVHolder(snvHolder);

				System.out.println("register " + flg + "\t" + flgtoStr(flg)
						+ "\t" + chr + "\t" + pos + "\t" + genomeR + "\t"
						+ normal.getRatio() + "\t" + tumor.getRatio());

	
			}
		}
		return ret;
	}

	private static String toTabstring(int[] sa) {

		StringBuffer sb = new StringBuffer();
		for (int n : sa) {
			sb.append(n + "\t");
		}
		return sb.toString();
	}

	private static String toTabstring(double[] sa) {

		StringBuffer sb = new StringBuffer();
		for (double d : sa) {
			sb.append(d + "\t");
		}
		return sb.toString();
	}

	private static int checkReg(PileUPResult normal, PileUPResult tumor,boolean targetSeq) {

		boolean regnormal = false;
		if (normal.isDiff()) {
			regnormal = normalCheck(normal);
		}
		boolean regTumor = false;
		if (tumor.isDiff()) {
			regTumor = tumorCheck(normal, tumor);
		}
		boolean depthcheckt = tumor.getTotal() >= KarkinosProp.mindepth;
		boolean depthcheckn = normal.getTotal() >= KarkinosProp.mindepthNormal;
		boolean depthcheck = depthcheckn && depthcheckt;

		if (!depthcheck) {

			// this block add 2012/10/18
			// need hetro snp info even where thre is
			// no tumor reads
			if (depthcheckn) {
				if (regnormal && regTumor) {
					return REGBOTH;
				} else if (regnormal) {
					return NormalSNP;
				}
			}
			return NONSignif;
		}

		if (regnormal && regTumor) {
			return REGBOTH;
		} else if (regnormal) {

			if (targetSeq) {

				if (normal.getRatio() < 0.35) {
					return NONSignif;
				}
			}
			return NormalSNP;
		} else if (regTumor) {
			return SomaticMutation;
		} else {
			return NONSignif;
		}

	}

	private static boolean tumorCheck(PileUPResult normal, PileUPResult tumor) {


		boolean tr = tumor.haveHigherRatio(KarkinosProp.min_initial_tumorratio);
		boolean nr = normal.haveLowerRatio(KarkinosProp.maxnormalratio);
		
		boolean rcheck = tr && nr;
		boolean morethanminsupport = true;
		if (!tumor.indel) {
			morethanminsupport = (tumor.getAltCnt() >= KarkinosProp.minsupportreads);
		} else {
			morethanminsupport = (tumor.indelcnt >= KarkinosProp.minsupportreads);
		}
		return rcheck && morethanminsupport;
	}

	private static boolean normalCheck(PileUPResult normal) {

		boolean rcheck = normal.getRatio() > KarkinosProp.normalSNPthres;
		// System.out.print("nr="+normal.getRatio());
		return rcheck;
	}

	final static int margin = 10;// for safety, theoretically not necessally

	private static int[] setPileUP(PileUPResult ret, String chr, int pos,
			char genomeR, List<SamHolder> samList, int startidx) {

		boolean debug = true;

		int[] reta = new int[3];
		boolean diff = false;
		ret.setGenomeRef(genomeR);
		IndelInfo indelinfo = new IndelInfo();
		int depth = 0;
		boolean isFirst = false;
		int retindex = startidx;
		int lowqual = 0;
		int ontarget = 0; 

		for (int n = startidx; n < samList.size(); n++) {

			SamHolder sh = samList.get(n);
			if(sh.getOi()!=null){
				if(sh.getOi().isOntag()){
					ontarget = 1;
				}
			}
			SAMRecord sam = sh.getSam();
			// if((pos ==
			// debugpos)&&sam.getReadName().equals("HWI-ST1181:113:D16B7ACXX:8:1301:21171:36898")){
			// System.out.println("here 2");
			// }
			if (!include(pos, sam)) {
				//
				if (sam.getAlignmentStart() > pos + margin) {
					// System.out.println("eidx"+aidx);
					break;
				}
				continue;
			}
			if (isFirst) {
				isFirst = false;
				retindex = n;
			}
			indelinfo.clear();
			int seqidx = getCharIdx(pos, sam, indelinfo);
			char ch = 0;
			byte qual = 0;
			int mapq = sam.getMappingQuality();
			if (seqidx >= 0) {

				try {
					ch = (char) sam.getReadBases()[seqidx];
					qual = sam.getBaseQualities()[seqidx];
				} catch (ArrayIndexOutOfBoundsException ex) {
					System.out.println(sam.format());
					ex.printStackTrace();
					continue;
				}

			}

			boolean diffread = isBase(ch) && ch != genomeR;

			if ((int) qual < KarkinosProp.minPhredQualForEach) {
				if (diffread) {
					lowqual++;
				}
				continue;
			}
			if (diffread) {
				diff = true;
			}
			ret.setBaseAndQual(ch, qual, mapq, indelinfo);
			depth++;
		}

		float lowqualratio = getRatio(lowqual, (lowqual + depth));
		ret.setLowqualratio(lowqualratio);

		ret.setDiff(diff);
		reta[0] = depth;
		reta[1] = retindex;
		reta[2] = ontarget;
		
		return reta;

	}

	private static float getRatio(int n, int m) {

		if (n == 0 || m == 0) {
			return 0f;
		}
		float r = (float) ((double) n / (double) m);
		return r;
	}

	private static boolean isBase(char ch) {
		return seqALL.indexOf(ch) >= 0;
	}

	private static int getCharIdx(int pos, SAMRecord sam, IndelInfo indelinfo) {

		indelinfo.indel = false;
		int start = sam.getAlignmentStart();
		int relpos = pos - start;

		int readidx = 0;
		int refidx = 0;

		if (relpos == 0) {
			List<CigarElement> list0 = sam.getCigar().getCigarElements();
			if (list0 != null && list0.size() >= 0) {
				CigarElement ce = list0.get(0);
				if (ce.getOperator() == (CigarOperator.SOFT_CLIP)) {
					return ce.getLength();
				}
			} else {
				return 0;
			}

		}

		List<CigarElement> list = sam.getCigar().getCigarElements();
		//
		if (list == null || list.size() == 0) {
			// Cigar does not exist which should not happen
			// assgin
			list = new ArrayList<CigarElement>();
			list.add(new CigarElement(sam.getReadLength(), CigarOperator.M));
		}

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
			
				if ((ce.getOperator() == CigarOperator.DELETION)||(ce.getOperator() == CigarOperator.N) ) {
					if (relpos == refidx + len) {
						indelinfo.indel = true;
						indelinfo.length = len;
						return -1;
					} else if (relpos < refidx + len) {
						return -1;
					}
				}
				if ((ce.getOperator() == CigarOperator.N) ) {
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

	private static String substring(String str, int s, int e) {

		if (e >= str.length()) {
			return str.substring(s);
		} else {
			return str.substring(s, e);
		}

	}

	private static boolean include(int pos, SAMRecord sam) {

		int start = sam.getAlignmentStart();
		int end = sam.getAlignmentEnd();
		if (end == 0 || end == start) {
			end = start + sam.getReadLength();
		}
		return pos >= start && pos <= end;
	}

}
