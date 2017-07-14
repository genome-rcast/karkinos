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
package jp.ac.utokyo.karkinos.cmdutils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class ExtractTP53GermLine {

	public static void main(String[] arg) throws IOException,
			NumberFormatException, MathException {

		// String indir =
		// "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast";

		String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/thcc,"
				+ "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/ncc,"
				+ "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/hcc-jp,"
				+ "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/goss," +
						"/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/tcga";

		List<File> snplist = new ArrayList<File>();
		for (String s : indir.split(",")) {
			File f = new File(s);
			searchRecursive2(f, snplist);
		}

		Map<Integer,SNPCounter> posmap = new HashMap<Integer,SNPCounter>();
		for (File f : snplist) {

			printSNP(f,posmap);

		}
//		Iterator<Integer> ite = posmap.keySet().iterator(); 
//		while(ite.hasNext()){
//			
//			int pos = ite.next();
//			System.out.println(pos+"\t"+posmap.get(pos).getId()+"\t"+posmap.get(pos).getCnt());
//			
//		}
		

	}

	private static void printSNP(File f, Map<Integer, SNPCounter> posmap) throws NumberFormatException,
			IOException, MathException {

		String id = f.getName().replaceAll("_normalsnp.vcf", "");
		int idx = id.indexOf(".vs");
		if (idx > 2) {
			id = id.substring(0, idx);
		}
		idx = id.indexOf("-THCC");
		if (idx > 2) {
			id = id.substring(0, idx);
		}
		System.out.print(id);
		String s = null;
		BufferedReader br = new BufferedReader(new FileReader(f));
		boolean print = false;
		String line = "";

		while ((s = br.readLine()) != null) {

			if (s.startsWith("#"))
				continue;
			String[] sa = s.split("\t");
			String info = sa[7];
			String[] infos = info.split(",");
			float af = getF(infos, "AF=");

			float bcn = getF(infos, "BCN=");
			boolean loh = false;
			loh = (bcn == 0f);

//			boolean hg17 = sa[0].equals("17");
//			int pos = Integer.parseInt(sa[1]);
//			boolean tp53 = (pos >= 7571720) && (pos <= 7590868);
//			boolean inr = (pos >= 7579450) && (pos <= 7579500);
			
			boolean hg17 = sa[0].equals("12");
			int pos = Integer.parseInt(sa[1]);
			boolean tp53 = (pos >= 112241766) && (pos <= 112241766);
			boolean inr = (pos >= 112241766) && (pos <= 112241766);
			
			
			boolean left = false;
			float aft = getF(infos, "AFT=");
			if (aft > 0.6 && (loh || af > 0.9)) {
				left = true;
			}

			String at = "GG";

			if (af >= 0.35) {

				at = "AG";
				if (af > 0.9) {
					at = "AA";
				}

				if (hg17 && tp53) {
					line = s;
					String str = " "+af+","+sa[0]+","+pos+","+1+","+sa[4]+"/"+sa[3]+",#"+sa[2];
					System.out.println(str);
					SNPCounter sq = null;
					if(!posmap.containsKey(pos)){
						sq = new SNPCounter();
						sq.setId(str);
						posmap.put(pos, sq);
					}else{
						sq = posmap.get(pos);
					}
					sq.setCnt(sq.getCnt()+1);
					
				}

				if (hg17 && inr) {

					print = true;
					//System.out.println(id + "\t" + at + "\t" + loh + "\t"
					//		+ left + "\t" + s);
				}

			}

		}

		if (print == false) {
			
				System.out.println();
			
			if (line != null) {
				String[] sa = line.split("\t");
				if (sa.length > 2) {
					String info = sa[7];
					String[] infos = info.split(",");

					float bcn = getF(infos, "BCN=");
					boolean loh = false;
					boolean left = false;
					loh = (bcn == 0f);

					//System.out.println(id + "\t" + "CC" + "\t" + loh + "\t"
					//		+ left + "\t" + line);
				}
			}
		}

	}

	private static float getF(String[] infos, String string) {

		float f = 0;
		for (String s : infos) {

			//
			if (s.startsWith(string)) {
				String ss = s.substring(s.indexOf("=") + 1);
				f = Float.parseFloat(ss);
			}

		}
		return f;
	}

	private static float exec(File f) throws IOException {

		float purity = 0f;
		String id = f.getName().replaceAll("_textdata.txt", "");

		String s = null;
		BufferedReader br = new BufferedReader(new FileReader(f));
		while ((s = br.readLine()) != null) {

			String[] sa = s.split("\t");
			if (s.startsWith("tc used")) {
				purity = Float.parseFloat(sa[1]);
			}

		}
		if (purity == 0.15)
			purity = 0.12f;
		if (id.equals("BN46T")) {
			purity = 0.701f;
		}
		return purity;

	}

	private static boolean searchRecursive2(File f, final List<File> list) {

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
					return searchRecursive2(new File(absPath), list);
				}
			}

		});

		if (listFiles == null)
			return false;
		for (File ff : listFiles) {
			if (ff.isFile() && ff.getName().contains("_normalsnp.vcf")) {

				// System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

}
