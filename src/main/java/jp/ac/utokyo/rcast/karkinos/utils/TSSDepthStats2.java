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
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class TSSDepthStats2 extends ReadWriteBase {

	public static void main(String[] arg) throws IOException {

		String bam = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/normal/"
				+ "Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_normal_genome.bam";
		
		//String bam = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/normal/"
		//	+ "Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_normal_genome.bam";

		String bandf = "/GLUSTER_DIST/data/users/ueda/project/Asada/ucscgene.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/project/Goto/RefSeqTSS_normal_500.txt";

		// /
		File f = new File("/GLUSTER_DIST/data/users/ueda/project/Goto/500.txt");

		SAMFileReader bamr = getReader(bam);

		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader sband = new BufferedReader(new FileReader(bandf));

		Map<String, String[]> mapucsc = new HashMap<String, String[]>();
		String line = null;
		while ((line = sband.readLine()) != null) {

			// if(cntn>10)break;
			if (line.startsWith("#")) {
				continue;
			}

			String[] sa = line.split("\t");

			String ucsc = sa[0];
			mapucsc.put(ucsc, sa);

		}

		// double[] forC = new double[2001];
		// double[] revC = new double[2001];

		Map<Integer, double[]> forCM = new HashMap<Integer, double[]>();
		Map<Integer, double[]> revCM = new HashMap<Integer, double[]>();
		for (int n = 0; n <= 6; n++) {

			forCM.put(n, new double[2001]);
			revCM.put(n, new double[2001]);

		}

		//
		BufferedReader two200 = new BufferedReader(new FileReader(f));
		line = null;
		int cnt[] = new int[7];
		while ((line = two200.readLine()) != null) {

			// if(cntn>10)break;
			if (line.startsWith("#")) {
				continue;
			}

			String[] sa = line.split("\t");
			int end = 6;
//			if(sa.length==2){
//				end = 1;
//			}
//			if(sa.length==3){
//				end = 2;
//			}
//			
			
			
			for (int n = 0; n <= end; n++) {

				//
				double[] forC = forCM.get(n);
				double[] revC = revCM.get(n);
	
				if(sa.length==1){
					forC = forCM.get(1);
					revC = revCM.get(1);
					boolean addf = addDepth(forC, revC, sa[1], mapucsc, bamr);
					if (addf) {
						cnt[1] = cnt[1]+1;
					}
					break;
				}
				
				boolean addf = addDepth(forC, revC, sa[n], mapucsc, bamr);
				if (addf) {
					cnt[n] = cnt[n]+1;
				}

			}

		}

		//

		for (int n = 0; n <= 6; n++) {

			System.out.println(cnt[n]);

		}

		for (int n = 0; n < 2000; n++) {

			for (int m = 0; m <= 6; m++) {

				double[] forC = forCM.get(m);
				double[] revC = revCM.get(m);

				bw.write(forC[n] + " \t" + revC[n] + " \t");

			}
			bw.write("\n");
		}
		bw.close();

	}

	private static boolean addDepth(double[] forC, double[] revC, String ucsc,
			Map<String, String[]> mapucsc, SAMFileReader bamr) {

		//
		if (!mapucsc.containsKey(ucsc)) {
			return false;
		}

		String[] sa = mapucsc.get(ucsc);
		if(sa.length<4){
			return false;
		}

		String chr = sa[1];
		int pos = Integer.parseInt(sa[3]);
		String stranbd = sa[2];
		int start = pos - 1000;
		int end = pos + 1000;

		boolean strand = sa[2].trim().equals("+");
		CloseableIterator<SAMRecord> ite = bamr.query(chr, start, end, false);
		List<ReadInterval> ri = getList(ite);

		if (strand) {

			for (int n = -1000; n <= 1000; n++) {
				int pos2 = pos + n;
				int[] ret = getD(pos2, ri);
				forC[n + 1000] = forC[n + 1000] + ret[0];
				revC[n + 1000] = revC[n + 1000] + ret[1];

			}

		} else {

			for (int n = 1000; n > -1000; n--) {
				int pos2 = pos + n;
				int[] ret = getD(pos2, ri);
				forC[n + 1000] = forC[n + 1000] + ret[1];
				revC[n + 1000] = revC[n + 1000] + ret[0];

			}

		}

		return true;

	}

	private static int[] getD(int pos2, List<ReadInterval> rin) {

		int f = 0;
		int r = 0;
		for (ReadInterval ri : rin) {

			if (ri.getStart() > pos2)
				continue;
			if (ri.getEnd() < pos2)
				continue;

			if (ri.isForward(pos2)) {
				f++;
			} else {
				r++;
			}

		}
		//
		return new int[] { f, r };

	}

	private static List<ReadInterval> getList(CloseableIterator<SAMRecord> itet) {

		Map<String, ReadInterval> map = new HashMap<String, ReadInterval>();

		while (itet.hasNext()) {

			SAMRecord sam = itet.next();
			if (sam.getReadUnmappedFlag())
				continue;
			if (!sam.getProperPairFlag())
				continue;
			//
			ReadInterval ri = null;
			if (map.containsKey(sam.getReadName())) {

				//
				ri = map.get(sam.getReadName());
				ri.add(sam);

			} else {

				ri = new ReadInterval(sam);
				map.put(sam.getReadName(), ri);

			}

		}
		itet.close();
		//
		//
		List<ReadInterval> list = new ArrayList<ReadInterval>();
		Iterator<String> ite = map.keySet().iterator();
		while (ite.hasNext()) {

			//
			list.add(map.get(ite.next()));

		}
		return list;

	}

}
