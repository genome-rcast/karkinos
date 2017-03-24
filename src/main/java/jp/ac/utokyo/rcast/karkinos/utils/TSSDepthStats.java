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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class TSSDepthStats extends ReadWriteBase {

	public static void main(String[] arg) throws IOException {

		String bam = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/normal/"
				+ "Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_normal_genome.bam";

		String bandf = "/GLUSTER_DIST/data/users/ueda/project/Asada/RefSeqTSS.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/project/Asada/RefSeqTSS_normal_ExpA.txt";

		// /
		File f = new File(
				"/data/users/ueda/project/Asada/OV_exp/IlluminaRNASeqV2_level3");
		Map<String, SummaryStatistics> map = new HashMap<String, SummaryStatistics>();

		int filecount = 0;
		for (File ff : f.listFiles()) {

			if (ff.getName().endsWith("isoforms.results")) {
				filecount++;
				// if(filecount>10) break;
				//
				BufferedReader br = new BufferedReader(new FileReader(ff));
				String line = br.readLine();// title
				while ((line = br.readLine()) != null) {

					//
					String[] sa = line.split("\t");
					String name = sa[0];
					Double val = 0d;
					try {
						val = Double.parseDouble(sa[2]);
					} catch (Exception ex) {

					}
					//
					SummaryStatistics ss = null;
					if (map.containsKey(name)) {
						//
						ss = map.get(name);
					} else {
						ss = new SummaryStatistics();
					}
					ss.addValue(val);
					map.put(name, ss);

				}
				br.close();
				System.out.println("filecount = " + filecount);
			}

		}

		//
		List<Double> valindex = new ArrayList<Double>();
		Iterator<String> ite0 = map.keySet().iterator();
		while (ite0.hasNext()) {
			//
			String s = ite0.next();
			double d = map.get(s).getGeometricMean();

			//
			if (d > 0) {
				valindex.add(d);
			}

		}
		Collections.sort(valindex);
		//
		int len = valindex.size() / 10;
		//
		double[] idxv = new double[10];
		for (int n = 0; n < 10; n++) {

			//
			idxv[n] = valindex.get(n * len);
			System.out.println(idxv[n]);

		}

		SAMFileReader bamr = getReader(bam);

		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader sband = new BufferedReader(new FileReader(bandf));

		// bw.write("pos/fr \t forward depth \t reverse depth \n");
		String line = null;
		//
		Set<String> kset = new HashSet<String>();

		// double[] forC = new double[2001];
		// double[] revC = new double[2001];

		Map<Integer, double[]> forCM = new HashMap<Integer, double[]>();
		Map<Integer, double[]> revCM = new HashMap<Integer, double[]>();
		//

		int cntn = 0;
		Map<Integer, Integer> counter = new HashMap<Integer, Integer>();

		while ((line = sband.readLine()) != null) {

			// if(cntn>10)break;
			if (line.startsWith("#")) {
				continue;
			}

			String[] sa = line.split("\t");

			String ucsc = sa[0];
			//
			if (!map.containsKey(ucsc)) {
				continue;
			}
			cntn++;

			// if(cntn>100)break;

			int idx = getIdx(map.get(ucsc), idxv);
			// System.out.println("index="+idx);
			if (counter.containsKey(idx)) {
				int n = counter.get(idx);
				n++;
				counter.put(idx, n);
			} else {
				counter.put(idx, 1);
			}

			double[] forC = null;
			double[] revC = null;

			if (forCM.containsKey(idx)) {

				//
				forC = forCM.get(idx);

			} else {
				forC = new double[2001];
				forCM.put(idx, forC);
			}
			if (revCM.containsKey(idx)) {

				//
				revC = revCM.get(idx);

			} else {
				revC = new double[2001];
				revCM.put(idx, revC);
			}

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

			// CloseableIterator<SAMRecord> ite = bamr.query(chr, start, end,
			// false);
			// List<SAMRecord> list = new ArrayList<SAMRecord>();
			// while (ite.hasNext()) {
			//
			// SAMRecord sam = ite.next();
			// list.add(sam);
			// }
			// ite.close();
			
			
			 CloseableIterator<SAMRecord> ite = bamr.query(chr, start, end,
			 false);
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

		}

		System.out.println(cntn);

		//
		System.out.println(counter);

		sband.close();
		bw.write("pos/fr \t");
		bw.write(" 0 forward depth \t 0 reverse depth \t");
		bw.write(" 1 forward depth \t 1 reverse depth \t");
		bw.write(" 2 forward depth \t 2 reverse depth \t");
		bw.write(" 3 forward depth \t 3 reverse depth \t");
		bw.write(" 4 forward depth \t 4 reverse depth \t");
		bw.write(" 5 forward depth \t 5 reverse depth \t");
		bw.write(" 6 forward depth \t 6 reverse depth \t");
		bw.write(" 7 forward depth \t 7 reverse depth \t");
		bw.write(" 8 forward depth \t 8 reverse depth \t");
		bw.write(" 9 forward depth \t 9 reverse depth \t");
		bw.write("\n");
		//
		for (int n = -1000; n <= 1000; n++) {

			bw.write(n + " \t");
			for (int m = 0; m <= 10; m++) {

				double[] forC = forCM.get(m);
				double[] revC = revCM.get(m);
				bw.write(forC[n + 1000] + " \t" + revC[n + 1000] + " \t");

			}
			bw.write(" \n");

		}

		bw.close();

	}

	private static int getIdx(SummaryStatistics ss, double[] idxv) {

		double d = ss.getGeometricMean();
		if (d == 0)
			return 0;
		int n = 0;
		for (; n < idxv.length; n++) {

			//
			double dd = idxv[n];
			if (dd > d)
				return n;

		}
		return n;
	}

	// private static int[] getD(int pos2, List<SAMRecord> list) {
	//
	// int f = 0;
	// int r = 0;
	// for (SAMRecord sam : list) {
	//
	// if (sam.getReadUnmappedFlag()) {
	// continue;
	// }
	// if(!sam.getFirstOfPairFlag()){
	// continue;
	// }
	// if (sam.getAlignmentStart() > pos2)
	// continue;
	// if (sam.getAlignmentEnd() < pos2)
	// continue;
	//
	// if (sam.getReadNegativeStrandFlag()) {
	//
	//
	// if(sam.getFirstOfPairFlag()){
	// r++;
	// }else{
	// f++;
	// }
	//
	// } else {
	//
	// if(sam.getFirstOfPairFlag()){
	// f++;
	// }else{
	// r++;
	// }
	//
	// }
	//
	// }
	// //
	// return new int[] { f, r };
	//
	// }

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
