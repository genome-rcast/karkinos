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
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.readssummary.CounterA;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.summary.CNAInterval;
import jp.ac.utokyo.rcast.karkinos.summary.ChromBand;

public class CheckHDAMPFromFile {

	public static void main(String[] arg) throws IOException, SQLException {

		//
		// String cnv = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0"
		// +
		// "/HCC_exome_rcast/THCC_12-7T-THCC_12-7N/THCC_12-7T-THCC_12-7N/THCC_12-7T-THCC_12-7N_cnvdepth.txt";
		// String cnvallelic =
		// "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0"
		// +
		// "/HCC_exome_rcast/THCC_12-7T-THCC_12-7N/THCC_12-7T-THCC_12-7N/THCC_12-7T-THCC_12-7N_cnvAllelicDepth.txt";

		String s4 = "/GLUSTER_DIST/data/users/ueda/genome/hg19chrband.txt";

		String ref = "/GLUSTER_DIST/data/Genomes/ucscgene_hg19/totoki_hg19/refFlat_120402.txt";
		GeneExons ex = new GeneExons(ref);

//		String out = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/karkinosExomeAmpHD.txt";
//
//		String out2 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/karkinosTotal.txt";
//
//		String out3 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/karkinosEachHigh.txt";
//
//		String out4 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/karkinosEachLow.txt";


		String out = "/GLUSTER_DIST/data/users/ueda/project/Mano/karkinosExomeAmpHD.txt";

		String out2 = "/GLUSTER_DIST/data/users/ueda/project/Mano/karkinosTotal.txt";

		String out3 = "/GLUSTER_DIST/data/users/ueda/project/Mano/karkinosEachHigh.txt";

		String out4 = "/GLUSTER_DIST/data/users/ueda/project/Mano/karkinosEachLow.txt";

		
		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);

		String census = "/GLUSTER_DIST/data/users/ueda/genome/cencercensus.txt";

		Set<String> cosmic = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(census));
		String geneline = null;
		while ((geneline = br.readLine()) != null) {
			String ge = geneline.split("\t")[0];
			cosmic.add(ge);
		}

		cosmic.add("SETDB1");
		cosmic.add("JAK3");
		cosmic.add("EYS");
		cosmic.add("MET");
		cosmic.add("VEGFA");

		br.close();
	//	String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/tcga";

//		 String indir =
//		 "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast,"
//		 +
//		 "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_NCC/karkinosver4.0.25_2/karkinosver4.0.25.130819,"
//		 +
//		 "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor/,"
//		 + "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/tcga,"
//		 + "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/goss,";

		// String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/test";
		
		 String indir = "/GLUSTER_DIST/data/users/ueda/project/Mano/karkinos_results";

		List<File> list = new ArrayList<File>();

		List<File> listtotal = new ArrayList<File>();

		ChromBand cb = new ChromBand(s4);

		// Map<String,Bean> map = new TreeMap<String,Bean>();
		Map<String, Bean> map = new LinkedHashMap<String, Bean>();

		Map<String, Map<String, Bean>> mapeach = new LinkedHashMap<String, Map<String, Bean>>();

		for (CNAInterval cnai : cb.getList()) {
			//
			// System.out.println(cnai.getName());
			map.put(cnai.getChr() + cnai.getName(), new Bean());
		}

		for (String s : indir.split(",")) {
			File f = new File(s);
			// chage to reccursive file search
			searchRecursive(f, list);
		}

		for (String s : indir.split(",")) {
			File f = new File(s);
			// chage to reccursive file search
			searchRecursive2(f, listtotal);
		}

		CheckHDAMPFromFile inst = new CheckHDAMPFromFile();

		List<String> hdList = new ArrayList<String>();
		List<String> ampList = new ArrayList<String>();
		Map<String, CounterA> mapgenehd = new HashMap<String, CounterA>();
		//
		Map<String, CounterA> mapgeneamp = new HashMap<String, CounterA>();

		for (File f : list) {

			File totalf = getTotalF(listtotal, f);
			try {
				inst.exec(f, pw, ex, map, cb, hdList, ampList, mapgenehd,
						mapgeneamp, totalf, mapeach);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		for (String sms : ampList) {

			String[] genes = sms.split("\t");
			String gene = "";
			if (genes.length > 8) {
				gene = genes[8];
			}
			String mostrec = getMostrec(gene, mapgeneamp, cosmic);
			pw.append(sms + "\t" + mostrec + "\n");
		}
		for (String hds : hdList) {

			String[] genes = hds.split("\t");
			String gene = "";
			if (genes.length > 8) {
				gene = genes[8];
			}
			String mostrec = getMostrec(gene, mapgenehd, cosmic);
			pw.append(hds + "\t" + mostrec + "\n");
		}

		pw.close();

		FileWriter fw2 = new FileWriter(new File(out2));
		BufferedWriter bw2 = new BufferedWriter(fw2);
		PrintWriter pw2 = new PrintWriter(bw2);
		Iterator<String> keys = map.keySet().iterator();

		while (keys.hasNext()) {

			String key = keys.next();
			Bean b = map.get(key);
			//
			if (b.high.getN() >= 1000) {
				pw2.print(key + "\t" + (float) b.high.getMean() + "\t"
						+ (float) b.low.getMean() + "\n");
			}
		}

		pw2.close();

		FileWriter fw3 = new FileWriter(new File(out3));
		BufferedWriter bw3 = new BufferedWriter(fw3);
		PrintWriter pw3 = new PrintWriter(bw3);
		Iterator<String> keys1 = mapeach.keySet().iterator();
		pw3.print("cb/sample" + "\t");
		while (keys1.hasNext()) {
			String key = keys1.next();
			pw3.print(key + "\t");
		}
		pw3.print("\n");
		//
		for (CNAInterval cnai : cb.getList()) {

			keys1 = mapeach.keySet().iterator();
			while (keys1.hasNext()) {
				String key = keys1.next();
				Map<String, Bean> me = mapeach.get(key);
				//
				String cbk = cnai.getChr() + cnai.getName();
				Bean b = me.get(cbk);
				if (b == null) {
					cbk = "chr"+cbk;
					b = me.get(cbk);
				}
				if (b == null) {
					cbk = cbk.replaceAll("chr", "");
					b = me.get(cbk);
				}
				if (b != null) {
					pw3.print(b.high.getMean() + "\t");
				} else {
					pw3.print(2 + "\t");
				}
			}
			pw3.print("\n");
		}

		pw3.close();

		FileWriter fw4 = new FileWriter(new File(out4));
		BufferedWriter bw4 = new BufferedWriter(fw4);
		PrintWriter pw4 = new PrintWriter(bw4);
		Iterator<String> keys2 = mapeach.keySet().iterator();
		pw4.print("cb/sample" + "\t");
		while (keys1.hasNext()) {
			String key = keys2.next();
			pw4.print(key + "\t");
		}
		pw4.print("\n");
		//
		for (CNAInterval cnai : cb.getList()) {

			keys2 = mapeach.keySet().iterator();
			while (keys2.hasNext()) {
				String key = keys2.next();
				Map<String, Bean> me = mapeach.get(key);
				//
				String cbk = cnai.getChr() + cnai.getName();
				Bean b = me.get(cbk);
				if (b == null) {
					cbk = "chr"+cbk;
					b = me.get(cbk);
				}
				if (b == null) {
					cbk = cbk.replaceAll("chr", "");
					b = me.get(cbk);
				}
				if (b != null) {
					pw4.print(b.low.getMean() + "\t");
				} else {
					pw3.print(1 + "\t");
				}
			}
			pw4.print("\n");
		}

		pw4.close();

	}

	private static File getTotalF(List<File> listtotal, File f) {

		String s = f.getName().replace("cnvAllelicDepth.txt", "");
		for (File ff : listtotal) {
			if (ff.getName().contains(s)) {
				return ff;
			}
		}
		return null;
	}

	private static String getMostrec(String genes, Map<String, CounterA> map,
			Set<String> cosmic) {

		//
		int max = 0;
		String maxs = "";
		for (String s : genes.split(",")) {

			if (s.length() > 2) {
				//
				if (cosmic.contains(s)) {
					if (map.containsKey(s)) {

						int cnt = map.get(s).getCnt();
						if (cnt > max) {
							max = cnt;
							maxs = s;
						}

					}
				}

			}

		}
		if (max > 1) {
			return maxs;
		}

		max = 0;
		maxs = "";
		for (String s : genes.split(",")) {

			if (s.length() > 2) {
				//
				if (map.containsKey(s)) {

					int cnt = map.get(s).getCnt();
					if (cnt > max) {
						max = cnt;
						maxs = s;
					}

				}

			}

		}
		if (max > 1) {
			return maxs;
		}

		return "";

	}

	private static boolean searchRecursive(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("cnvAllelicDepth.txt");
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
			if (ff.isFile() && ff.getName().contains("cnvAllelicDepth.txt")) {

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
						|| name.contains("cnvdepth.txt");
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
			if (ff.isFile() && ff.getName().contains("cnvdepth.txt")) {

				System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

	public void exec(File f, PrintWriter pw, GeneExons ex,
			Map<String, Bean> map, ChromBand cb, List<String> hdList2,
			List<String> ampList2, Map<String, CounterA> mapgenehd,
			Map<String, CounterA> mapgeneamp, File totalf,
			Map<String, Map<String, Bean>> mapeach) throws IOException,
			SQLException {

		List<CopyNumberInterval> copyNumberIntervalListTotal = getTotal(totalf);

		String id = f.getName();
		id = id.replaceAll("_cnvAllelicDepth.txt", "");

		Map<String, Bean> map2 = new HashMap<String, Bean>();
		mapeach.put(id, map2);

		// cnv
		BufferedReader br = new BufferedReader(new FileReader(f));
		String s = null;
		int cnt = 0;
		List<CopyNumberInterval> copyNumberIntervalList = new ArrayList<CopyNumberInterval>();

		float aafb4 = 0f;
		float bafb4 = 0f;
		String chrb4 = "";
		CopyNumberInterval cni = null;
		SummaryStatistics ss = new SummaryStatistics();

		while ((s = br.readLine()) != null) {

			if (s.startsWith("#"))
				continue;
			try {
				cnt++;
				String[] sa = s.split("\t");
				String chr0 = sa[0];
				if (!chr0.contains("chr")) {
					chr0 = "chr" + chr0;
				}
				float aaf = Float.parseFloat(sa[9]);
				float baf = Float.parseFloat(sa[10]);

				if (aaf == 3) {
					ss.addValue((double) (Math.abs(Float.parseFloat(sa[5]) - 1)));
				}
				if (aaf == 4) {
					ss.addValue((double) (Math.abs(Float.parseFloat(sa[5]) - 1) / 2));
				}
				if (baf == 1) {
					ss.addValue((double) (Math.abs(1 - Float.parseFloat(sa[6]))));
				}

				int pos = Integer.parseInt(sa[1]);

				float highs = Float.parseFloat(sa[5]);
				float lows = Float.parseFloat(sa[6]);
				if (highs > 4)
					highs = 4;
				if (lows > 4)
					lows = 4;
				// System.out.println(highs+"\t"+lows);

				//
				// String key = chr0+"-"+pos;
				String key = chr0 + cb.getBand(chr0, pos, pos);
				//
				Bean b = null;
				if (map.containsKey(key)) {
					b = map.get(key);
				} else {
					b = new Bean();
					map.put(key, b);
				}
				b.high.addValue(highs);
				b.low.addValue(lows);

				Bean b2 = null;
				if (map2.containsKey(key)) {
					b2 = map2.get(key);
				} else {
					b2 = new Bean();
					map2.put(key, b);
				}
				b2.high.addValue(highs);
				b2.low.addValue(lows);

				//
				boolean statesame = (chrb4.equals(chr0)) && (aaf == aafb4)
						&& (baf == bafb4);
				boolean statechange = !statesame;

				aafb4 = aaf;
				bafb4 = baf;
				chrb4 = chr0;
				//
				if (statechange) {

					//
					if (cni != null) {
						copyNumberIntervalList.add(cni);
					}
					cni = new CopyNumberInterval();
					cni.setChr(chr0);
					cni.setStart(pos);
					cni.setEnd(pos);
					cni.setAaf(aaf);
					cni.setBaf(baf);
				} else {

					if (cni != null) {
						cni.setEnd(pos);
					}
				}
			} catch (Exception ex2) {
				continue;
			}

		}
		if (cni != null) {
			copyNumberIntervalList.add(cni);
		}
		List<CopyNumberInterval> ampList = new ArrayList<CopyNumberInterval>();
		List<CopyNumberInterval> hdList = new ArrayList<CopyNumberInterval>();

		for (CopyNumberInterval cni0 : copyNumberIntervalList) {

			//
			int len = cni0.getEnd() - cni0.getStart();
			if (len < 5000000) {

				if (cni0.getAaf() >= 3 && cni0.getBaf() >= 2) {
					// System.out.println(cni0.getChr() + "\t" + cni0.getAaf()
					// + "\t" + cni0.getBaf()+"\t"+len);
					//
					ampList.add(cni0);
				}
				if (cni0.getAaf() <= 1 && cni0.getBaf() <= 1) {

					hdList.add(cni0);

				}

			}
		}

		List<CopyNumberInterval> ampListCheck = new ArrayList<CopyNumberInterval>();
		double peakdist = ss.getMean();
		if (ss.getN() < 10 || peakdist < 0.5) {

			peakdist = 0.5;

		}
		//
		for (CopyNumberInterval cni1 : ampList) {

			//
			double checkok = checkOK(cni1, f, peakdist, ampList);
			if (checkok > 2) {

				cni1.setCopynumber((float) checkok);
				ampListCheck.add(cni1);

			}

		}

		for (CopyNumberInterval cni1 : ampListCheck) {

			String chr = cni1.getChr();
			if (!chr.contains("chr")) {
				chr = "chr" + chr;
			}
			int len = cni1.getEnd() - cni1.getStart();
			String gss = ex.getGeneSymbols(chr, cni1.getStart(), cni1.getEnd());

			String band = chr
					+ cb.getBand(chr, cni1.getStart(), cni1.getStart());

			String amps = "AMP" + "\t" + band + "\t" + id + "\t" + chr + "\t"
					+ cni1.getStart() + "\t" + cni1.getEnd() + "\t"
					+ cni1.getCopynumber() + "\t" + len + "\t" + gss;
			System.out.println(amps);

			setGene(mapgeneamp, gss);

			ampList2.add(amps);

			// pw.print("AMP" + "\t" + id + "\t" + chr + "\t"
			// + cni1.getStart() + "\t" + cni1.getEnd() + "\t"
			// + cni1.getCopynumber() + "\t" + len +"\t"+gss+ "\n");

		}
		List<CopyNumberInterval> hdListCheck = new ArrayList<CopyNumberInterval>();
		for (CopyNumberInterval cni2 : hdList) {

			int cnthd = count(cni2.getChr(), hdList);
			//check chr count
			if (cnthd < 4) {

				 if(checkCNTotal(cni2,copyNumberIntervalListTotal)){
				
	
					hdListCheck.add(cni2);

				 }
			}

		}

		for (CopyNumberInterval cni2 : hdListCheck) {

			String chr = cni2.getChr();
			if (!chr.contains("chr")) {
				chr = "chr" + chr;
			}

			int len = cni2.getEnd() - cni2.getStart();
			String gss = ex.getGeneSymbols(chr, cni2.getStart(), cni2.getEnd());
			//

			String band = chr
					+ cb.getBand(chr, cni2.getStart(), cni2.getStart());

			String hds = "HD" + "\t" + band + "\t" + id + "\t" + chr + "\t"
					+ cni2.getStart() + "\t" + cni2.getEnd() + "\t"
					+ cni2.getCopynumber() + "\t" + len + "\t" + gss;
			System.out.println(hds);
			hdList2.add(hds);

			setGene(mapgenehd, gss);

			// pw.print("HD" + "\t" + id + "\t" + chr + "\t"
			// + cni2.getStart() + "\t" + cni2.getEnd() + "\t"
			// + cni2.getCopynumber() + "\t" + len +"\t"+gss+ "\n");

		}

	}

	private boolean checkCNTotal(CopyNumberInterval cni,
			List<CopyNumberInterval> copyNumberIntervalListTotal) {

		float cn = getMinCN(cni, copyNumberIntervalListTotal);
		
		int maxcn =1;
		if(cn4>cn2){
			maxcn=2;
		}
		return cn <= maxcn;
	}

	private float getMinCN(CopyNumberInterval cni,
			List<CopyNumberInterval> copyNumberIntervalListTotal) {

		if (copyNumberIntervalListTotal.size() == 0) {
			return 0;
		}

		float mincn = 2.0f;
		for (CopyNumberInterval cn : copyNumberIntervalListTotal) {

			//
			if (cni.getChr().equals(cn.getChr())) {

				if ((cni.getStart() <= cn.getEnd())
						&& (cn.getStart() <= cni.getEnd())) {
					// contain
					if (cn.getCopynumber() < mincn) {
						//
						mincn = cn.getCopynumber();
					}

				}

			}

		}
		return mincn;
	}
	
	int cn4=0;
	int cn2=0;
	private List<CopyNumberInterval> getTotal(File totalf) throws IOException {

		cn4=0;
		cn2=0;
		
		List<CopyNumberInterval> copyNumberIntervalList = new ArrayList<CopyNumberInterval>();
		if (totalf == null) {
			return copyNumberIntervalList;
		}

		BufferedReader br = new BufferedReader(new FileReader(totalf));
		String s = null;
		int cnt = 0;
		CopyNumberInterval cni = null;
		String chrb4 = "";
		float cnb4 = 0;
		float aafb4 = 0;
		float bafb4 = 0;

		while ((s = br.readLine()) != null) {

			if (s.startsWith("#"))
				continue;
			try {
				cnt++;
				String[] sa = s.split("\t");
				String chr0 = sa[0];
				if (!chr0.contains("chr")) {
					chr0 = "chr" + chr0;
				}
				int pos = Integer.parseInt(sa[1]);
				float cn = Float.parseFloat(sa[10]);
				float aaf = Float.parseFloat(sa[8]);
				float baf = Float.parseFloat(sa[9]);

				if(cn==2.0f){
					cn2++;
				}else if(cn==4.0f){
					cn4++;
				}
				
				//
				boolean statesame = (chrb4.equals(chr0)) && (aaf == aafb4)
						&& (baf == bafb4);

				boolean statechange = !statesame;

				cnb4 = cn;
				chrb4 = chr0;
				aafb4 = aaf;
				bafb4 = baf;

				//
				if (statechange) {

					//
					if (cni != null) {
						copyNumberIntervalList.add(cni);
					}
					cni = new CopyNumberInterval();
					cni.setChr(chr0);
					cni.setStart(pos);
					cni.setEnd(pos);
					cni.setCopynumber(cn);
					cni.setAaf(aaf);
					cni.setBaf(baf);

				} else {

					if (cni != null) {
						cni.setEnd(pos);
					}
				}
			} catch (Exception ex2) {
				continue;
			}

		}
		if (cni != null) {
			copyNumberIntervalList.add(cni);
		}
		return copyNumberIntervalList;

	}

	private void setGene(Map<String, CounterA> map, String gss) {

		for (String s : gss.split(",")) {

			if (s.length() > 2) {

				CounterA c = null;
				if (map.containsKey(s)) {
					c = map.get(s);
				} else {
					c = new CounterA();
					map.put(s, c);
				}
				c.inc();
			}

		}

	}

	private double checkOK(CopyNumberInterval cni1, File f, double peakdist,
			List<CopyNumberInterval> ampList) throws IOException {

		int countchr = count(cni1.getChr(), ampList);
		if (countchr > 4) {
			return -1;
		}

		BufferedReader br = new BufferedReader(new FileReader(f));
		String s = null;
		while ((s = br.readLine()) != null) {

			if (s.startsWith("#"))
				continue;
			try {
				String[] sa = s.split("\t");
				String chr0 = sa[0];
				if (!chr0.contains("chr")) {
					chr0 = "chr" + chr0;
				}
				int pos = Integer.parseInt(sa[1]);
				float high = Float.parseFloat(sa[7]);
				//
				if (cni1.getChr().equals(chr0)) {

					//
					if ((pos >= cni1.getStart()) && (pos <= cni1.getEnd())) {
						//
						if (high > 1 + (peakdist * 2)) {

							return (high / peakdist) + 2;

						}

					}

				}
			} catch (Exception ex) {
				continue;
			}
		}
		return -1;

	}

	private int count(String chr, List<CopyNumberInterval> ampList) {
		int n = 0;

		for (CopyNumberInterval cni : ampList) {

			if (chr.equals(cni.getChr())) {
				n++;
			}

		}

		return n;
	}

}
