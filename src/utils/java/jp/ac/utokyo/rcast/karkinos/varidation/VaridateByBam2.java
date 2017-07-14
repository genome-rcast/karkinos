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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPPool;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class VaridateByBam2 extends ReadWriteBase {

	public static void main(String[] arg) throws Exception {

		
		String s ="/data/users/ueda/project/Aihara/20140516Editome/editingpos.txt";
		String out ="/data/users/ueda/project/Aihara/20140516Editome/resultout.txt";
		
		BufferedReader snvsb = new BufferedReader(new FileReader(s));
		String line = null;
		//
		List<String> posS = new ArrayList<String>();
		while ((line = snvsb.readLine())!= null) {
			
			//
			posS.add(line);	
			
		}
		
		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		
		
		//
		File f = new File("/data/users/yamamoto/RNAseq/rnaseq/Glioma/bam");
		
		Map<String,Map<String,PileUPResult>> mm
			= new HashMap<String,Map<String,PileUPResult>>();
		
		List<String> fnl = new ArrayList<String>(); 
		
		for(File ff:f.listFiles()){
			
			if(ff.getName().endsWith(".bam")){
				
				System.out.println(ff.getName());
				fnl.add(ff.getName());
				
				Map<String,PileUPResult> m = new HashMap<String,PileUPResult>();
				mm.put(ff.getName(),m);
				
				SAMFileReader bamr = getReader(ff);
				for(String ss:posS){
					
					String[] sa = ss.split("\t");
					String chr = sa[0];
					int pos = Integer.parseInt(sa[1]);
					//
					char nuc = tgr.getGenomeNuc(chr, pos, true);
					String ref =  Character.toUpperCase(nuc)+"";
					String alt = "G";
					if(ref.equals("T")){
						alt="C";
					}
					
					PileUPResult pur = getData(bamr, chr,pos,
							ref,alt);
					m.put(ss, pur);
					
				}
				bamr.close();
				//break;
			}
			
		}
		
		//print
		bw.append("cand/name \t");
		for(String ns:fnl){
			
			bw.append(ns);
			bw.append("\t");
			
		}
		bw.append("\n");
		
		for(String ss:posS){
			
			bw.append(ss);
			bw.append("\t");
			for(String ns:fnl){
				
				PileUPResult pur = mm.get(ns).get(ss);
				bw.append(toAF(pur)+"");
				bw.append("\t");
				
			}
			bw.append("\n");
			
		}
		
		bw.close();
		
	}

	private static float toAF(PileUPResult pur) {
		// TODO Auto-generated method stub
		float r =  pur.getRatio();
		if(r<0.01){
			r= 0;
		}
		return r;
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

	private static PileUPResult getData(SAMFileReader bamr, String chr, int pos,
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

	private static int setPileup(String chr, int pos, char genomeR, PileUPResult ret,
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
