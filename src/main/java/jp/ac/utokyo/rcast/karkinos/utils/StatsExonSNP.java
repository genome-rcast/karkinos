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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.summary.Filebean;

public class StatsExonSNP {

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

		//String indir = "/GLUSTER_DIST/data/users/yamamoto/exome";
		String indir = 
			"/GLUSTER_DIST/data/users/yamamoto/exome," +	
			"/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast," +
				"/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_ncc";
		//String outpath = "/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP3.vcf";
		String outpath = "/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP5.vcf";

//		String indir ="/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor";
//		String outpath="/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNPvc.vcf";
//		
		List<File> list = new ArrayList<File>();
		Set<String> nameset = new HashSet<String>();

		for (String s : indir.split(",")) {
			File f = new File(s);
			// chage to reccursive file search
			searchRecursive(f, list, nameset);
		}

		Map<String, Map<Integer, SNPbean>> map = new LinkedHashMap<String, Map<Integer, SNPbean>>();

		//
		int samplesize = list.size();
		for (File f : list) {

			//
			if(f.getName().contains("THCC_135T")){
				continue;
			}
			if(f.getName().contains("THCC_12-15T")){
				continue;
			}
			if(f.getName().contains("THCC_150T")){
				continue;
			}
			addData(f, map);

		}
		output(map, samplesize, outpath);

	}

	private static void addData(File f, Map<String, Map<Integer, SNPbean>> map)
			throws IOException {

		//
		FileInputStream fis = new FileInputStream(f);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));

		try {

			int totalcnt = 0;
			int init = fis.available();
			String chr = "";
			for (;;) {
				int point = init - fis.available();
				String line = br.readLine();
				if (line == null)
					break;
				if (line.startsWith("#")) {
					continue;
				}
				totalcnt++;
				String[] sa = line.split("\t");
				String _chr = sa[0];
				chr = _chr;
				int pos = Integer.parseInt(sa[1]);
				String id = sa[2];

				// non dbSNP exon SNP
				// register
				Map<Integer, SNPbean> chrData = map.get(chr);
				if (chrData == null) {
					chrData = new TreeMap<Integer, SNPbean>();
					map.put(chr, chrData);
				}
				//
				SNPbean snpBean = chrData.get(pos);
				if (snpBean == null) {

					snpBean = new SNPbean(sa);
					chrData.put(pos, snpBean);

				} else {
					snpBean.inc(sa);
				}

			}

		} finally {
			br.close();
		}

	}

	private static void output(Map<String, Map<Integer, SNPbean>> map,
			int samplesize, String outpath) {

		//
		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			// header
			bw.write("#chr" + "\t");
			bw.write("pos" + "\t");
			bw.write("ID" + "\t");
			bw.write("ref" + "\t");
			bw.write("alt" + "\t");
			bw.write("count" + "\t");
			bw.write("ratio" + "\t");
			bw.write("meanAF" );
			bw.write("\n");

			writeCNV(bw, map, samplesize);

			bw.close();
		} catch (IOException ex) {

		}

	}

	private static void writeCNV(BufferedWriter bw,
			Map<String, Map<Integer, SNPbean>> map, int total)
			throws IOException {

		Iterator<String> ite = map.keySet().iterator();
		while (ite.hasNext()) {

			String chr = ite.next();
			Map<Integer, SNPbean> snpData = map.get(chr);
			Iterator<Integer> ite2 = snpData.keySet().iterator();
			while (ite2.hasNext()) {
				
				int pos = ite2.next();
				SNPbean bean = snpData.get(pos);
				bw.write(chr + "\t");
				bw.write(pos + "\t");
				bw.write(bean.id + "\t");
				bw.write(bean.ref + "\t");
				bw.write(bean.alt + "\t");
				bw.write(bean.cnt + "\t");
				bw.write(ratio(bean.cnt,total) + "\t");
				bw.write(bean.getMean()+"");
				bw.write("\n");
			}

		}

	}

	private static float ratio(int cnt, int total) {
		
		return (float)((double)((double)cnt/(double)total));
	}

	private static boolean searchRecursive(File f, final List<File> list,
			final Set<String> nameset) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("_normalsnp.vcf");
				if (!target)
					return false;
				//
				String absPath = dir.getAbsolutePath() + File.separator + name;
				if (new File(absPath).isFile()) {
					return true;
				} else {
					return searchRecursive(new File(absPath), list, nameset);
				}
			}

		});
		if (listFiles == null)
			return false;

		List<File> flist = new ArrayList<File>();
		for (File ff : listFiles) {
			if (ff.isFile()) {
				flist.add(ff);
			}
		}

		for (File ff : flist) {

			if (ff.getName().contains("_normalsnp.vcf")) {

				if (nameset.contains(ff.getName())) {

				} else {
					list.add(ff);
					nameset.add(ff.getName());
				}

			}
		}
		return true;

	}
}
