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
package jp.ac.utokyo.rcast.karkinos.varidation;

import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPPool;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;

public class FilterByBam extends ReadWriteBase {

	public static void main(String[] arg) {

		String s1 = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/challenge4/challenge4submit8.vcf";

		String bam1 = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/synthetic.challenge.set3.normal_m.bam";
		String bam2 = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/"
				+ "normal/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_normal_genome.bam";
		String bam3 = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/synthetic.challenge.set1.normal.v2.bam";
		String bai = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/synthetic.challenge.set1.normal.v2.bai";
		
		String s2 = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/challenge4/challenge4submit10.vcf";

		try {
			SAMFileReader bamr3 = getReader(bam3);
			BAMIndexer indexer = new BAMIndexer(new File(bai), bamr3.getFileHeader());

			bamr3.enableFileSource(true);
			int totalRecords = 0;

			// create and write the content
			for (SAMRecord rec : bamr3) {
				System.out.println(rec.format());
				if (++totalRecords % 1000000 == 0) {
					System.out.println(totalRecords + " reads processed ...");
				}
				indexer.processAlignment(rec);
			}
			System.out.print(totalRecords);
			indexer.finish();
			bamr3.close();
			//FilterByBam._main(s1, s2, bam1, bam2, bam3);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void _main(String s1, String s2, String bam1, String bam2,
			String bam3) throws IOException {

		FilterByBam inst = new FilterByBam();
		inst.exec(s1, s2, bam1, bam2, bam3);

	}

	//
	private void exec(String snvs, String out, String bam1, String bam2,
			String bam3) throws IOException {

		//
		BufferedWriter br = new BufferedWriter(new FileWriter(out));
		BufferedReader snvsb = new BufferedReader(new FileReader(snvs));
		String line = snvsb.readLine();// title
		br.write(line);

		//
		while (line != null) {

			line = snvsb.readLine();
			if (line == null)
				break;
			if (line.startsWith("#")) {

				br.write(line + "\n");
				continue;
			}
			//
			String[] data = line.split("\t");

			String chr = data[0];
			int pos = 0;

			try {
				pos = Integer.parseInt(data[1]);
			} catch (Exception ex) {
				continue;
			}

			// if(pos==79914583){
			// System.out.println("stop");
			// }else{
			// continue;
			// }
			String ref = data[3];
			String alt = data[4];

			String qual = data[5];
			int qi = Integer.parseInt(qual);

			//
			if (isNumber(chr)) {
				chr = "chr" + chr;
			}
			if (chr.contains("X") || chr.contains("Y") || chr.contains("M")) {
				chr = "chr" + chr;
			}

			System.out.println(line);

			SAMFileReader bamr1 = getReader(bam1);

			PileUPResult pr = getData(bamr1, chr, pos, ref, alt);

			int thres = 3;
			float thresf = 0.04f;

			boolean twosup = (qi < 210);
			if (twosup) {
				thres = 2;
				thresf = 0.025f;
			}
			twosup = (qi < 120);
			if (twosup) {
				thres = 1;
				thresf = 0.015f;
			}

			//
			boolean erroevalidated = (pr.getRefAltCnt()[1] >= thres)
					&& (pr.getRatio() > thresf);
			if (erroevalidated) {
				System.out.println("by1");
				continue;
			}

			SAMFileReader bamr2 = getReader(bam2);

			pr = getData(bamr2, chr, pos, ref, alt);
			erroevalidated = (pr.getRefAltCnt()[1] >= thres)
					&& (pr.getRatio() > thresf);
			if (erroevalidated) {
				System.out.println("by2");
				continue;
			}

			SAMFileReader bamr3 = getReader(bam3);
			pr = getData(bamr3, chr, pos, ref, alt);
			erroevalidated = (pr.getRefAltCnt()[1] >= 2)
					&& (pr.getRatio() > 0.03);
			if (erroevalidated) {
				System.out.println("by3");
				continue;
			}

			br.write(line + "\n");
		}
		br.close();

	}

	private boolean isNumber(String chr) {

		try {
			Integer.parseInt(chr);
		} catch (Exception ex) {
			return false;
		}
		return true;
	}

	private PileUPResult getData(SAMFileReader bamr, String chr, int pos,
			String ref, String alt) {

		//
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		SAMRecordIterator ite = bamr.query(chr, pos, pos, false);
		while (ite.hasNext()) {

			//
			SAMRecord sam = ite.next();
			reads.add(sam);

		}
		ite.close();
		boolean isindel = (ref.length() > 1) || (alt.length() > 1);

		PileUPResult pr = PileUPPool.borrowObject();
		char genomeR = ref.charAt(0);
		pr.setGenomeRef(genomeR);
		setPileup(chr, pos, genomeR, pr, reads, isindel);
		return pr;

	}

	private int setPileup(String chr, int pos, char genomeR, PileUPResult ret,
			List<SAMRecord> reads, boolean isindel) {

		int maxreadsend = -1;
		boolean diff = false;

		int depth = 0;

		for (SAMRecord sam : reads) {

			// System.out.println(sam.format());
			IndelInfo indelinfo = new IndelInfo();
			int seqidx = getCharIdx(pos, sam, indelinfo);
			char ch = 0;
			byte qual = 0;
			int readslen = sam.getReadLength();
			int readsend = 0;

			if ((seqidx >= 0) && (seqidx < readslen)) {
				ch = (char) sam.getReadBases()[seqidx];
				qual = sam.getBaseQualities()[seqidx];
				readsend = Math.min(readslen - seqidx, seqidx);
				if (indelinfo.indel) {
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

	private static String substring(String str, int s, int e) {

		if (e >= str.length()) {
			return str.substring(s);
		} else {
			return str.substring(s, e);
		}

	}

}
