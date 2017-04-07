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
package jp.ac.utokyo.rcast.karkinos.utils;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class BInStats extends ReadWriteBase {

	public static void main(String[] arg) throws IOException {

		String bam = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/normal/"
				+ "Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_normal_genome.bam";
		
		
		// String bandf =
		// "/GLUSTER_DIST/data/users/ueda/genome/hg19chrband.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/project/Asada/Se68-cfDNA_25ng-Se68-bDNA_25ng_normal_genome.1m.txt";

		SAMFileReader bamr = getReader(bam);

		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		// BufferedReader sband = new BufferedReader(new FileReader(bandf));

		bw.write("chr \t start \t end \t name \t gieStain \t forward \t reverse \n");
		String line = null;

		//
		SAMSequenceDictionary ssd = bamr.getFileHeader()
				.getSequenceDictionary();
		for (SAMSequenceRecord ssr : ssd.getSequences()) {

			String chr = ssr.getSequenceName();
			int len = ssr.getSequenceLength();
			// String chr = sa[0];
			// int start = Integer.parseInt(sa[1]);
			// int end = Integer.parseInt(sa[2]);
			//
			int mod = len / 1000000;
			for (int m = 0; m <= mod; m++) {

				int start = m*1000000;
				int end = (m+1)*1000000;
				if(end>len){
					end = len;
				}
				
				int forward = 0;
				int reverse = 0;
				CloseableIterator<SAMRecord> ite = bamr.query(chr, start, end,
						true);
				while (ite.hasNext()) {

					SAMRecord sam = ite.next();
					if (sam.getReadUnmappedFlag()) {
						continue;
					}
					if (sam.getReadNegativeStrandFlag()) {
						reverse++;
					} else {
						forward++;
					}

				}
				ite.close();
				String outstr = chr + "\t" + start + "\t" + end + "\t"
						+ forward + "\t"
						+ reverse;
				bw.write(outstr + "\n");
				System.out.println(outstr);
			}
		}
		bw.close();

	}

}
