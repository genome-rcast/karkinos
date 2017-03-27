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
import java.util.List;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class ExtractGermLine {

	public static void main(String[] arg) throws IOException, NumberFormatException, MathException {

	//	String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast";
		
		String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/thcc,"
			+ "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/ncc,"
			+ "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/hcc-jp,"
			//
			+ "/GLUSTER_DIST/data/users/ueda/ICGC/karkinosResult/goss";
		
		String out = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/germlineSNP.txt";
		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);

		//
		List<File> list = new ArrayList<File>();
		for (String s : indir.split(",")) {
			File f = new File(s);
			searchRecursive(f, list);
		}

		Map<String, Float> map = new HashMap<String, Float>();
		for (File f : list) {

			String id = f.getName().replaceAll("_textdata.txt", "");
			float purity = 0.12f;
			try {
				purity = exec(f);
			} catch (IOException e) {
				e.printStackTrace();
			}
			// System.out.println(id+"\t"+purity);
			map.put(id, purity);
		}
		//
		//

		List<File> snplist = new ArrayList<File>();
		for (String s : indir.split(",")) {
			File f = new File(s);
			searchRecursive2(f, snplist);
		}

		for (File f : snplist) {

			String id = f.getName().replaceAll("_normalsnp.vcf", "");
			float purity = 0.5f;
			if(map.containsKey(id)){
				purity = map.get(id);
			}
			// System.out.println(id+"\t"+purity);
			writeSNP(f, pw, purity, id);

		}
		pw.close();
	}

	private static void writeSNP(File f, PrintWriter pw, float purity, String id)
			throws NumberFormatException, IOException, MathException {

		String s = null;
		BufferedReader br = new BufferedReader(new FileReader(f));
		SummaryStatistics ss = new SummaryStatistics();
		
		while ((s = br.readLine()) != null) {

			if (s.startsWith("#"))
				continue;

			String[] sa = s.split("\t");
			String info = sa[7];
			String[] infos = info.split(",");
			float aaf = getF(infos, "ACN=");
			float baf = getF(infos, "BCN=");
			float depth = getF(infos, "DP=");
			//
			if (aaf == baf) {

				if (depth >= 20) {

					float af = getF(infos, "AF=");
					if (af >= 0.35 && af <= 0.65) {

						float aft = getF(infos, "AFT=");
						//System.out.println(s);
						double ratio = aft/af;
						//System.out.println(af + "\t" + aft +"\t"+ratio);
						ss.addValue(ratio);
					}

				}

			}

		}
		
		double mean = ss.getMean();
		double sd = ss.getStandardDeviation();
		
		NormalDistribution nd = new NormalDistributionImpl(mean,sd);
		double alpha = 0.05; 
		double upper = nd.inverseCumulativeProbability(1-alpha/2);
		double lower = nd.inverseCumulativeProbability(alpha/2);
		
		br = new BufferedReader(new FileReader(f));
		while ((s = br.readLine()) != null) {
			
			
			if (s.startsWith("#"))
				continue;
			String[] sa = s.split("\t");
			String info = sa[7];
			String[] infos = info.split(",");
			float aaf = getF(infos, "ACN=");
			float baf = getF(infos, "BCN=");
			float depth = getF(infos, "DP=");
			//
			if (baf==0) {

				if (depth >= 20) {

					float af = getF(infos, "AF=");
					if (af >= 0.35 && af <= 0.65) {

						float aft = getF(infos, "AFT=");
						//System.out.println(s);
						double ratio = aft/af;
						if(ratio>upper){
							pw.println(id+"\t"+af+"\t"+aft+"\t"+aaf+"\t"+baf+"\t"+s);
							System.out.println(id+"\t"+af+"\t"+aft+"\t"+aaf+"\t"+baf+"\t"+s);
						}
						if(ratio<lower){
							pw.println(id+"\t"+af+"\t"+aft+"\t"+aaf+"\t"+baf+"\t"+s);
							System.out.println(id+"\t"+af+"\t"+aft+"\t"+aaf+"\t"+baf+"\t"+s);
						}
						
						
					}

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

	private static boolean searchRecursive(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("_textdata.txt");
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
			if (ff.isFile() && ff.getName().contains("_textdata.txt")) {

				// System.out.println(ff.getName());
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
