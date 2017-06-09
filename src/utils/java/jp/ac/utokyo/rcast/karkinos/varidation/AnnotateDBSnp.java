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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPPool;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;

public class AnnotateDBSnp extends ReadWriteBase {

	public static void main(String[] arg) {

			try {
				AnnotateDBSnp._main(arg);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

	}

	

	public static void _main(String[] arg) throws IOException {

		
		AnnotateDBSnp inst = new AnnotateDBSnp();
		inst.exec();
		
	}

	//
	private void exec() throws IOException {

		//
		String snvs = "/GLUSTER_DIST/data/users/ueda/project/sakata/data.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/project/sakata/dataadd3.txt";
		
		String dbSNP = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_snp132.txt";
		String g1000 = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_ALL.sites.2012_02.txt";
		String cosmic  = "/usr/local/karkinos/karkinos/genome/CosmicCompleteExport_v62_291112.vcf";
		String exome = "/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP5.vcf";
		
		DbSNPAnnotation dba = new DbSNPAnnotation(null, null,0.01f,null,exome);
		
		BufferedWriter br = new BufferedWriter(new FileWriter(out));
		BufferedReader snvsb = new BufferedReader(new FileReader(snvs));
		String line = snvsb.readLine();
	
		
		//
		while (line != null) {

			line = snvsb.readLine();
			if(line==null)break;
			line = line.trim();
			System.out.println("line "+ line);
			
			
			if (line == null)
				break;
			if (line.startsWith("#")) {

				br.write(line + "\n");
				continue;
			}
			//
			String[] data = line.split("\t");
			if(line.length()==0){
				br.write( "\n");
				continue;
			}
			String chr = data[0];
			int pos = 0;
			
			try{
				pos = Integer.parseInt(data[1]);
			}catch(Exception ex){
				br.write( "\n");
				continue;
			}
			DbSNPBean dbs = dba.lookup(chr, pos);

			
			String s = "";
			if(dbs!=null){
				s= dbs.getInfo();				
			}
			
			String freq = "";
			if(dbs!=null){
				freq= dbs.getFreq()+"";				
			}

						
			br.write(line + "\t" +s+"\t"+freq + "\n");


		}
		br.close();
		snvsb.close();

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
