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

public class StatErrorRate {

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

//		// String indir = "/GLUSTER_DIST/data/users/yamamoto/exome";
//		String indir = "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.10/MT";
//		// String outpath = "/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP3.vcf";
//		// String outpath = "/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP4.vcf";
//
//		// String indir
//		// ="/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor";
//		// String outpath="/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNPvc.vcf";
//		//
//		List<File> list = new ArrayList<File>();
//		Set<String> nameset = new HashSet<String>();
//
//		for (String s : indir.split(",")) {
//			File f = new File(s);
//			// chage to reccursive file search
//			searchRecursive(f, list, nameset);
//		}
//
//		//
//		int samplesize = list.size();
//		for (File f : list) {
//
//			int[] ret = getRet(f);
//			System.out.println(f.getName() +"\t"+ ret[0] + "\t" + ret[1] +"\t"+ (double)ret[0]/(double)(ret[0]+ret[1]));
//
//		}
		//
		StatErrorRate inst = new StatErrorRate();
		inst.exec();	
		
	}
	
	public void exec() throws IOException{
		
		String f = "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.10" +
		"/MT/MT5_1-MT5N/MT5_1-MT5N_annover_input.txt.genome_summaryannoplus.csv";

//
		Map<String,Bean> map = getRET(new File(f));
		//
		Iterator<String> ite = map.keySet().iterator();
		while(ite.hasNext()){
			
			String key = ite.next();
			Bean b = map.get(key); 
			//System.out.println("aaa");
			String s = key+"\t"+b.getRatio();
			System.out.println(s);
			
		}
	
	}
	
	class Bean{
		
		int pass;
		int error;
		float getRatio(){
			if(error==0){
				if(pass==0)return 0f;
				return 1;
			}
			return (float)((double)pass/(double)(pass+error));
		}
		
	}
	
	private Map<String,Bean> getRET(File f) throws IOException {

		//
		Map<String,Bean> map = new HashMap<String,Bean>();
		FileInputStream fis = new FileInputStream(f);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		int pass = 0;
		int error = 0;
		try {
			 br.readLine();//title
			String l = null;
			while ((l = br.readLine()) != null) {

				String[] sa = l.split(",");
				if (sa.length > 10) {
					
					int depth = Integer.parseInt(sa[23]);
					float ratio = Float.parseFloat(sa[25]);
					String key = getKey(ratio,depth);
					Bean b = null;
					if(map.containsKey(key)){
						b = map.get(key);
					}else{
						b = new Bean();
						map.put(key, b);
					}
					
					if (sa[8].equals("PASS")) {
						b.pass = b.pass+1;				
						
					} else {
						b.error = b.error+1;	
					}
				}

			}

		} finally {
			br.close();
		}		
		return map;
		
	}


	private static String getKey(float ratio, int depth) {
		
		String s1 = String.valueOf(((int)(ratio * 10)));
		String s2 = String.valueOf(((int)(Math.log(depth)/Math.log(2))));
		return s1+"-"+s2;
		
	}

	private static int[] getRet(File f) throws IOException {

		//
		FileInputStream fis = new FileInputStream(f);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		int pass = 0;
		int error = 0;
		try {

			String l = null;
			while ((l = br.readLine()) != null) {

				String[] sa = l.split(",");
				if (sa.length > 10) {
					if (sa[8].equals("PASS")) {
						pass++;
					} else {
						error++;
					}
				}

			}

		} finally {
			br.close();
		}
		return new int[] { pass, error };

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
			bw.write("meanAF");
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
				bw.write(ratio(bean.cnt, total) + "\t");
				bw.write(bean.getMean() + "");
				bw.write("\n");
			}

		}

	}

	private static float ratio(int cnt, int total) {

		return (float) ((double) ((double) cnt / (double) total));
	}

	private static boolean searchRecursive(File f, final List<File> list,
			final Set<String> nameset) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("annoplus.csv");
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

			if (ff.getName().contains("annoplus.csv")) {

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
