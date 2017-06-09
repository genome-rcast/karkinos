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
package jp.ac.utokyo.rcast.karkinos.summary.cnvbygene;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class CNVByGene {

	public static void main(String[] arg) throws IOException {

		String gr = "/usr/local/karkinos/karkinos/genome/refgenes_anno.csv";
		Map<String, GeneBean> map = CNVUtils.getMap(gr);

		String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast";
		List<File> depthlist = new ArrayList<File>();
		for (String s : indir.split(",")) {
			File f = new File(s);
			searchRecursive(f, depthlist);
		}

		List<File> depthAlist = new ArrayList<File>();
		for (String s : indir.split(",")) {
			File f = new File(s);
			searchRecursive2(f, depthAlist);
		}
		//

		for (File dfile : depthlist) {

			File allelicf = getAllelicF(depthAlist, dfile);
			Iterator<String> ite = map.keySet().iterator();
			Map<String, TreeMap<Integer, CNVEachBean>> depthmap = getDepthMap(
					dfile, true);

			Map<String, TreeMap<Integer, CNVEachBean>> adepthmap = getDepthMap(
					allelicf, false);

			float peakdistT[] = getdist(depthmap,true);
			float peakdistA[] = getdist(adepthmap,false);
			
			while (ite.hasNext()) {

				GeneBean b = map.get(ite.next());
				// System.out.println(b.getName() + "\t" + b.getChr() + "\t"
				// + b.getStart() + "\t" + b.getEnd());

			}
		}

	}

	
	
	private static float[] getdist(
			Map<String, TreeMap<Integer, CNVEachBean>> map, boolean total) {
		
		//
		float[] fa = new float[2];
		
		Iterator<String> ite = map.keySet().iterator();
		int cn2cnt=0;
		int cn4cnt=0;
		
		SummaryStatistics ss = new SummaryStatistics();
		
		while(ite.hasNext()){
			
			String chr = ite.next();
			TreeMap<Integer, CNVEachBean> tp = map.get(chr);
			Iterator<Integer> posite = tp.keySet().iterator();			
			
			while(posite.hasNext()){
				
				int pos = posite.next();
				CNVEachBean bean = tp.get(pos);
				if(total){
					
					float ar = bean.getRatioAdj();
					float cn = bean.getCopyNumber();
					if(cn==4.0f){
						cn4cnt++;
					}else if(cn==2.0f){
						cn2cnt++;
					}
					double diff = 0;
					if(cn!=2){
						
						if(cn>2){
							diff = (ar-1)/(cn-2);
						}else{
							diff = (ar-1)/(cn-2);
						}
						
					}							
					
					
					
					
				}else{
					
					
					
					
					
				}
				int baseploidy = cn2cnt>cn4cnt?2:4;
				//fa[0]= ss.getMean();
				//fa[1]= baseploidy;
			}
			
		}
		

		return fa;
	}



	private static Map<String, TreeMap<Integer, CNVEachBean>> getDepthMap(File f,
			boolean total) throws IOException {

		//
		Map<String, TreeMap<Integer, CNVEachBean>> map = new HashMap<String, TreeMap<Integer, CNVEachBean>>();
		//
		BufferedReader br = new BufferedReader(new FileReader(f));
		String s = null;
		while ((s = br.readLine()) != null) {

			//
			if (s.startsWith("#"))
				continue;
			try {

				String[] sa = s.split("\t");
				String chr0 = sa[0];
				if (!chr0.contains("chr")) {
					chr0 = "chr" + chr0;
				}
				int pos = Integer.parseInt(sa[1]);
				TreeMap<Integer, CNVEachBean> data = null;
				//
				if(map.containsValue(chr0)){
					data = map.get(chr0);
				}else{
					data = new TreeMap<Integer, CNVEachBean>();
					map.put(chr0, data);
				}
		
				CNVEachBean bean = new CNVEachBean(sa,total);
				data.put(pos, bean);
			
			} catch (Exception ex) {

			}
		}
		return map;

	}

	private static File getAllelicF(List<File> listtotal, File f) {

		String s = f.getName().replace("_cnvdepth.txt", "");
		for (File ff : listtotal) {
			if (ff.getName().contains(s)) {
				return ff;
			}
		}
		return null;
	}

	private static boolean searchRecursive(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("_cnvdepth.txt");
				if (!target)
					return false;
				//
				String absPath = dir.getAbsolutePath() + File.separator + name;
				if (new File(absPath).isFile()) {
					return true;
				} else {
					return searchRecursive(new File(absPath), list);
				}
			}

		});

		if (listFiles == null)
			return false;
		for (File ff : listFiles) {
			if (ff.isFile() && ff.getName().contains("_cnvdepth.txt")) {

				System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

	private static boolean searchRecursive2(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("_cnvAllelicDepth.txt");
				if (!target)
					return false;
				//
				String absPath = dir.getAbsolutePath() + File.separator + name;
				if (new File(absPath).isFile()) {
					return true;
				} else {
					return searchRecursive2(new File(absPath), list);
				}
			}

		});

		if (listFiles == null)
			return false;
		for (File ff : listFiles) {
			if (ff.isFile() && ff.getName().contains("_cnvAllelicDepth.txt")) {

				System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

}
