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
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class VaridateByBam extends ReadWriteBase {

	public static void main(String[] arg) {

		// String[] sa = new String[] {
		// "189","195","197","301","305","314","318",
		// "319","320","324","327","332","338","344","356","357","360","368"
		// ,"372","379","381","385","387","389","392","401","404","408",
		// "410","412","415","417","418","424","425","427","431","437","438"};

		// String[] sa = new String[] {
		// "7","27","35","49","50","54","58","65","94","131","167",
		// "172","199","200","317","325","326","329","337","351","352",
		// "362","362","378","380","394","405","424","426","434","445",
		// "475","507","511","512","516","521","524","530","531","532",
		// "533","534","537","538","548","552","553","562","566","567",
		// "570","571","576","581","582"
		// };

		// String[] sa = new String[] { "183"};
		//
		String[] sa = new String[] { "567" };

		for (String s : sa) {

			// String s1 =
			// "/GLUSTER_DIST/data/users/yamamoto/work/midorikawa/EHCC_karkinos/karkinos4.1.8/"
			// +"M_E_"+s+"T-M_E_"+s+"N/"
			// +
			// "M_E_"+s+"T-M_E_"+s+"N_annover_input.txt.genome_summaryannoplus.csv";

			// String s1 =
			// "/GLUSTER_DIST/data/users/yamamoto/work/midorikawa/EHCC_karkinos/karkinos4.1.8/"
			// +"M_E_"+s+"_T-M_E_"+s+"_N/"
			// +
			// "M_E_"+s+"_T-M_E_"+s+"_N_annover_input.txt.genome_summaryannoplus.csv";

			// String s1
			// ="/GLUSTER_DIST/data/users/yamamoto/work/midorikawa/EHCC_karkinos/karkinos4.1.8/M_E_567T_2-M_E_567N/"
			// +
			// "M_E_567T_2-M_E_567N_annover_input.txt.genome_summaryannoplus.csv"
			// ;
			//
			// String s2 =
			// "/GLUSTER_DIST/data/users/yamamoto/RNAseq/rnaseq/MHCC/MHCC_RNAseq/M_stR_"
			// + s
			// + "T-2/"
			// + "Aligned/M_stR_"
			// + s
			// + "T-2_100/M_stR_"
			// + s
			// + "T-2_100_fine.bam";
			//
			// String s3 = "/GLUSTER_DIST/data/users/ueda/pgtest/ehcc/";

			String s1 = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.11/Se-68_tuDNA-Se-68_b-DNA/" +
					"Se-68_tuDNA-Se-68_b-DNA_annover_input.txt.genome_summaryannoplus.csv";
			
//			String s2 = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.8/summary_cfDNA1/bam/"
//					+ "Se-67-cf25-Se-67_b-DNA_tumor_genome.bam";

//			String s2 = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se-71_ascites-Se-71_b-DNA/tumor/"
//				+"Se-71_ascites-Se-71_b-DNA_tumor_genome.bam";
			
			String s2 = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/tumor/"
				+ "Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_tumor_genome.bam";
		
			
			// String s1 =
			// "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.8/Se-68-tumor25-Se-68_b-DNA/"
			// +"Se-68-tumor25-Se-68_b-DNA_annover_input.txt.genome_summaryannoplus.csv";
			//
			// String s2 =
			// "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.8/summary_cfDNA1/bam/"+
			// "Se-68-cf10-Se-68_b-DNA_tumor_genome.bam";

			String s3 = "/GLUSTER_DIST/data/users/ueda/project/Asada/";

			System.out.println(s1);

			List<String> l = new ArrayList<String>();
			add(l, "-snvs", s1);
			add(l, "-bam", s2);
			add(l, "-out", s3);

			String[] ar = l.toArray(new String[l.size()]);

			File f1 = new File(s1);
			if (!f1.exists()) {
				continue;
			}
			System.out.println(s2);

			File f2 = new File(s2);
			if (!f2.exists()) {
				continue;
			}

			try {
				VaridateByBam._main(ar);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

	private static List<Option> getOption() {
		List<Option> optionlist = new ArrayList<Option>();

		optionlist.add(getOption("snvs", "snvList", true,
				"list Of SNV Candidates", true));

		optionlist.add(getOption("bam", "bam", true,
				"bam file used for varidation", true));

		optionlist.add(getOption("out", "out", true, "output directory", true));

		return optionlist;
	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	public static void _main(String[] arg) throws IOException {

		BasicParser parcer = new BasicParser();
		List<Option> optionList = getOption();
		Options opts = new Options();
		for (Option opt : optionList) {
			opts.addOption(opt);
		}

		CommandLine cl = null;
		try {
			cl = parcer.parse(opts, arg);
		} catch (ParseException e1) {
			System.out.println(e1.getMessage());
			HelpFormatter help = new HelpFormatter();
			help.setOptionComparator(new OptionComparator(optionList));
			help.printHelp("karkinos.jar varidate", opts, true);
			return;
		}
		//
		String snvs = null;
		if (cl.hasOption("snvs")) {
			snvs = cl.getOptionValue("snvs");
		}
		//
		String bam = null;
		if (cl.hasOption("bam")) {
			bam = cl.getOptionValue("bam");
		}
		//
		String out = null;
		if (cl.hasOption("out")) {
			out = cl.getOptionValue("out");
		}

		//
		VaridateByBam inst = new VaridateByBam();
		inst.exec(snvs, bam, out);

	}

	//
	private void exec(String snvs, String bam, String out) throws IOException {

		//
		out = out + (new File(snvs)).getName() + "validate.csv";
		BufferedWriter br = new BufferedWriter(new FileWriter(out));
		BufferedReader snvsb = new BufferedReader(new FileReader(snvs));
		String line = snvsb.readLine();
		br.write(line
				+ "QUAL2,TReads,QUAL3,vREF,vALT,vRatio,vTOTAL,Varidated \n");
		SAMFileReader bamr = getReader(bam);
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
			String[] data = line.split(",");
			String chr = data[1];
			int pos = 0;
			
			try{
				pos = Integer.parseInt(data[2]);
			}catch(Exception ex){
				continue;
			}
			
//			if(pos==79914583){
//				System.out.println("stop");
//			}else{
//				continue;
//			}
			String ref = data[4];
			String alt = data[5];

			boolean indel = (ref.length() > 1 || alt.length() > 1);
			if (indel) {
				System.out.println(line);
			}
			//
			if (isNumber(chr)) {
				chr = "chr" + chr;
			}
			if (chr.contains("X")||chr.contains("Y")||chr.contains("M")) {
				chr = "chr" + chr;
			}
			
			
			PileUPResult pr = getData(bamr, chr, pos, ref, alt);
			pr.setRatiocalclated(false);
			if (indel) {
				
				pr.setIndel(true);
			}else{
				pr.setIndel(false);
			}
			//
			boolean validated = (pr.getRefAltCnt()[1] >= 2)
					&& (pr.getRatio() > 0.05);
			if (indel) {
				validated = (pr.getRefAltCnt()[1] >= 1)
						&& (pr.getRatio() > 0.01);
			}
			
			if (pr.getRatio() > 0.3) {
				validated = true;
			}
			br.write(line + pr.getRefAltCnt()[0] + "," + pr.getRefAltCnt()[1]
					+ "," + pr.getRatio() + "," + pr.getTotalcnt() + ","
					+ validated + "\n");
//			if(ref.length()>1){
//			System.out.println( pr.getRefAltCnt()[0] + "," + pr.getRefAltCnt()[1]
//					+ "," + pr.getRatio() + "," + pr.getTotalcnt() + ","
//					+ validated + "\n");
	//		}
			//

		}
		bamr.close();
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

			//System.out.println(sam.format());
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
