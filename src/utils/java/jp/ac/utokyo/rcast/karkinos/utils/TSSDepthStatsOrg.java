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
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class TSSDepthStatsOrg extends ReadWriteBase {

	public static void main(String[] arg) throws IOException {

		String bam = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68-cfDNA_25ng-Se68-bDNA_25ng/normal/"
				+ "Se68-cfDNA_25ng-Se68-bDNA_25ng_normal_genome.bam";

		String bandf = "/GLUSTER_DIST/data/users/ueda/project/Asada/RefSeqTSS.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/project/Asada/RefSeqTSS_normal_TERT.txt";

		// /
		

		SAMFileReader bamr = getReader(bam);

		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader sband = new BufferedReader(new FileReader(bandf));

		// bw.write("pos/fr \t forward depth \t reverse depth \n");
		String line = null;
		//
		Set<String> kset = new HashSet<String>();

		 double[] forC = new double[2001];
		 double[] revC = new double[2001];

	
		//

		int cntn = 0;
		while ((line = sband.readLine()) != null) {

			// if(cntn>100)break;
			if (line.startsWith("#")) {
				continue;
			}

			String[] sa = line.split("\t");

			String ucsc = sa[0];
			//
			
			if(!ucsc.equals("uc003jcb.1")){
				continue;
			}
		
			cntn++;
			
			//if(cntn>100)break;
		

			String chr = sa[1];
			int pos = Integer.parseInt(sa[5]);
			String stranbd = sa[3];
			int start = pos - 1000;
			int end = pos + 1000;

			String s = sa[7];
			// if(s.equals("NR_046018")){
			// System.out.println("here");
			// }else{
			// continue;
			// }

			boolean strand = sa[2].trim().equals("+");

			String key = chr + "_" + pos;

			if (kset.contains(key)) {
				continue;
			}
			kset.add(key);

			System.out.println(line);

			//
			int forward = 0;
			int reverse = 0;
			CloseableIterator<SAMRecord> ite = bamr.query(chr, start, end,
					false);
			List<SAMRecord> list = new ArrayList<SAMRecord>();
			while (ite.hasNext()) {

				SAMRecord sam = ite.next();
				list.add(sam);
			}
			ite.close();

			if (strand) {

				for (int n = -1000; n <= 1000; n++) {
					int pos2 = pos + n;
					int[] ret = getD(pos2, list);
					forC[n + 1000] = forC[n + 1000] + ret[0];
					revC[n + 1000] = revC[n + 1000] + ret[1];

				}

			} else {

				for (int n = 1000; n > -1000; n--) {
					int pos2 = pos + n;
					int[] ret = getD(pos2, list);
					forC[n + 1000] = forC[n + 1000] + ret[1];
					revC[n + 1000] = revC[n + 1000] + ret[0];

				}

			}
			// System.out.println(outstr);
		}
		sband.close();
		bw.write("pos/fr \t");
		bw.write(" 0 forward depth \t 0 reverse depth \t");
		bw.write("\n");
		//
		for (int n = -1000; n <= 1000; n++) {

			bw.write(n+" \t");

			
				bw.write(forC[n+1000]+" \t"+revC[n+1000]+" \t");
				

			bw.write(" \n");

		}

		bw.close();

	}

	private static int getIdx(SummaryStatistics ss, double[] idxv) {

		double d = ss.getGeometricMean();
		if(d==0)return 0;
		int n = 0;
		for (; n < idxv.length; n++) {

			//
			double dd = idxv[n];
			if (dd > d)
				return n;

		}
		return n;
	}

	private static int[] getD(int pos2, List<SAMRecord> list) {

		int f = 0;
		int r = 0;
		for (SAMRecord sam : list) {

			if (sam.getReadUnmappedFlag()) {
				continue;
			}
			if (sam.getAlignmentStart() > pos2)
				continue;
			if (sam.getAlignmentEnd() < pos2)
				continue;

			if (sam.getReadNegativeStrandFlag()) {
				r++;
			} else {
				f++;
			}

		}
		//
		return new int[] { f, r };

	}

}
