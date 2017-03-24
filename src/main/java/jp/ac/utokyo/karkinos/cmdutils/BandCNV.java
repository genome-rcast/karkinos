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
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.summary.CNAInterval;
import jp.ac.utokyo.rcast.karkinos.summary.ChromBand;

public class BandCNV {

	public static void _main(String[] arg) throws IOException, SQLException {

		//
		String cbhg19f = "/GLUSTER_DIST/data/users/ueda/genome/hg19chrband.txt";
		String cbhg18f = "/GLUSTER_DIST/data/users/ueda/genome/hg18chrband.txt";
		ChromBand cbhg18 = new ChromBand(cbhg18f, "hg18");
		ChromBand cbhg19 = new ChromBand(cbhg19f, "hg19");
		// String out = "/GLUSTER_DIST/data/users/ueda/baylorCNV.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/baylorCNVsummary.txt";
		//

		String array = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/SNP6/selected/Baylor/txt";
		String exome = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast";

		// List<File> list = new ArrayList<File>();
		// searchRecursive(new File(array), list);

		Map<String, Map<String, Bean>> mmap = new LinkedHashMap<String, Map<String, Bean>>();

		for (File f : new File(array).listFiles()) {

			//
			if (f.getName().contains("allele.cnt")) {
				System.out.println(f.getName());

				String id = f.getName().replaceAll("_allele.cnt", "");
				id = id.replaceAll("PJ_Liver", "HCC-JP");

				Map<String, Bean> map = exec(id, f, cbhg18);
				mmap.put(id, map);
			}

		}
		//

		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);

		Iterator<String> ite = mmap.keySet().iterator();
		pw.print("" + "\t");
		while (ite.hasNext()) {
			pw.print(ite.next() + "\t");
		}
		pw.print("\n");

		for (CNAInterval cna : cbhg18.getList()) {

			String id0 = cna.getChr() + cna.getName();
			Iterator<String> ite2 = mmap.keySet().iterator();
			pw.print(id0 + "\t");
			while (ite2.hasNext()) {

				String id = ite2.next();
				Map<String, Bean> mapArray = mmap.get(id);
				Bean arrayB = mapArray.get(id0);
				float highA = 1;
				float lowA = 1;
				if (arrayB != null && arrayB.high.getN() > 100) {
					highA = (float) arrayB.high.getMean();
					lowA = (float) arrayB.low.getMean();
				}
				int flg = getFlg(highA, lowA);
				pw.print(flg + "\t");
			}
			pw.print("\n");

		}

		pw.close();

	}

	public static void main(String[] arg) throws Exception {

		//
		String cbhg19f = "/GLUSTER_DIST/data/users/ueda/genome/hg19chrband.txt";
		// String cbhg18f =
		// "/GLUSTER_DIST/data/users/ueda/genome/hg18chrband.txt";
		// ChromBand cbhg18 = new ChromBand(cbhg18f, "hg18");
		ChromBand cbhg19 = new ChromBand(cbhg19f, "hg19");
		// String out = "/GLUSTER_DIST/data/users/ueda/baylorCNV.txt";
		// String out = "/GLUSTER_DIST/data/users/ueda/earlyCNVsummary.txt";
		// String out =
		// "/GLUSTER_DIST/data/users/ueda/project/Aihara/AIsummary.txt";
		//String out = "/GLUSTER_DIST/data/users/ueda/project/Mano/PCNSL_summary.txt";
		//
		String out = "/GLUSTER_DIST/data/users/ueda/earlyCNVsummaryEHCC.txt";

		// String array =
		// "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/SNP6/selected/Baylor/txt";
		// String exome =
		// "/GLUSTER_DIST/data/users/yamamoto/exome/EHCC/EHCC_exome";
		 String exome =
		 "/GLUSTER_DIST/data/users/ueda/project/Midorikawa/EHCC";

		//String exome = "/GLUSTER_DIST/data/users/ueda/project/Mano/karkinos_results/";
		
		//String exome = "/GLUSTER_DIST/data/users/ueda/project/Mano/karkinos_results";
		
		List<File> list = new ArrayList<File>();
		searchRecursive(new File(exome), list);
		
		//
		

		List<File> list2 = new ArrayList<File>();
		searchRecursive2(new File(exome), list2);
		Map<String, Float> tp = getTP(list2);

		Map<String, List<String>> hdAmp = getHdAmp(list2);

		Map<String, Map<String, Bean>> mmap = new LinkedHashMap<String, Map<String, Bean>>();

		for (File f : list) {

			//

			String id = f.getName().replaceAll("_cnvAllelicDepth", "");

			try {
				Map<String, Bean> map = exec2(id, f, cbhg19, tp);
				mmap.put(id, map);
			} catch (Exception ex) {
			}

		}
		//

		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);

		Iterator<String> ite = mmap.keySet().iterator();
		pw.print("" + "\t");
		while (ite.hasNext()) {
			pw.print(ite.next() + "\t");
		}

		pw.print("lowave \t");
		pw.print("highave \t");

		pw.print("\n");

		for (CNAInterval cna : cbhg19.getList()) {

			String id0 = cna.getChr() + cna.getName();
			Iterator<String> ite2 = mmap.keySet().iterator();
			pw.print(id0 + "\t");

			SummaryStatistics hiave = new SummaryStatistics();
			SummaryStatistics lowave = new SummaryStatistics();
			//

			while (ite2.hasNext()) {

				String id = ite2.next();
				Map<String, Bean> mapArray = mmap.get(id);
				Bean arrayB = mapArray.get(id0);
				float highA = 1;
				float lowA = 1;
				if (arrayB != null && arrayB.high.getN() > 30) {
					highA = (float) arrayB.high.getMean();
					lowA = (float) arrayB.low.getMean();
					hiave.addValue(highA);
					lowave.addValue(lowA);

				}
				List<String> hdAm =  hdAmp.get(id);
				int flg = getFlg2(highA, lowA,hdAm,cna);
				pw.print(flg + "\t");
				// pw.print(lowA+":"+highA+"\t");
			}

			pw.print(lowave.getMean() + "\t");
			pw.print(hiave.getMean());

			pw.print("\n");

		}

		pw.close();

	}

	private static Map<String, List<String>> getHdAmp(List<File> list2) throws Exception {

		Map<String, List<String>> m = new HashMap<String, List<String>>();

		for (File f : list2) {
			
			List<String> l = new ArrayList<String>();
			BufferedReader abr = new BufferedReader(new FileReader(f));
			float tp = 0f;
			String as = null;
			while ((as = abr.readLine()) != null) {

				if(as.endsWith("HD")){
					l.add(as);
				}
				if(as.endsWith("AMP")){
					l.add(as);
				}

			}

			if (tp < 0.4)
				tp = 0.4f;

			// System.out.println(f.getName());
			String id = f.getName().replaceAll("_textdata", "");

			m.put(id, l);

		}
		return m;

	}

	private static Map<String, Float> getTP(List<File> list2) throws Exception {

		Map<String, Float> m = new HashMap<String, Float>();

		for (File f : list2) {
			BufferedReader abr = new BufferedReader(new FileReader(f));
			float tp = 0f;
			String as = null;
			while ((as = abr.readLine()) != null) {

				if (as.startsWith("tc used")) {
					tp = Float.parseFloat(as.split("\t")[1]);
				}

			}

			if (tp < 0.4)
				tp = 0.4f;

			// System.out.println(f.getName());
			String id = f.getName().replaceAll("_textdata", "");

			m.put(id, tp);

		}
		return m;
	}

	private static int getFlg(float highA, float lowA) {

		if (highA == 1 && lowA == 1) {
			return -1;
		}

		if (highA < 0.4 && lowA > -0.4) {
			return 0;
		}

		if (highA >= 0.4 && lowA <= -0.4) {
			return 3;
		}

		if (highA >= 0.4) {
			return 1;
		}

		if (lowA <= -0.4) {
			return 2;
		}

		return 0;
	}

	private static int getFlg2(float highA, float lowA, List<String> hdAm, CNAInterval cna) {

		if (highA == 1 && lowA == 1) {
			return -1;
		}

		if (highA < 1.75 && lowA > 0.75) {
			return 0;
		}

		if (highA >= 1.75 && lowA <= 0.75) {
			return 9;
		}

		if (highA >= 3.25) {
			if (containFocal(hdAm, cna,true)) {
				return 4;
			}
			return 3;
		}
		if (highA >= 2.5) {
			return 2;
		}
		if (highA >= 1.75) {
			return 1;
		}

		if (lowA <= 0.75) {

			if (lowA <= 0.5) {
				if (containFocal(hdAm, cna,false)) {
					return 6;
				}
			}
			return 7;
		}

		return 0;
	}

	private static boolean containFocal(List<String> hdAm,
			CNAInterval cna,boolean gain) {
		//
		for(String s:hdAm){
			//
			String[] sa = s.split("\t");
			String chr = sa[0];
			int start = Integer.parseInt(sa[1].replaceAll(",",""));
			int end = Integer.parseInt(sa[2].replaceAll(",",""));
			
			//
			if(s.endsWith("HD")&&gain){
				continue;
			}
			if(s.endsWith("AMP")&&!gain){
				continue;
			}
			//
			if(intercect(cna,chr,start,end)){
				return true;
			}
			
			
		}
		return false;
	}

	private static boolean intercect(CNAInterval cna, String chr, int start,
			int end) {
		
		//
		boolean bool1 = cna.getChr().equals(chr);
		
		if(!bool1) return false;
		//				
		return cna.intercect(start, end);
	
	}

	private static Map<String, Bean> exec(String id, File arrayf,
			ChromBand cbhg18) throws IOException, SQLException {

		Map<String, Bean> mapArray = new LinkedHashMap<String, Bean>();

		BufferedReader abr = new BufferedReader(new FileReader(arrayf));
		String as = null;
		while ((as = abr.readLine()) != null) {
			String[] asa = as.split("\t");
			// float lows = Float.parseFloat(asa[5]);
			// float highs = Float.parseFloat(asa[6]);

			float lows = Float.parseFloat(asa[5]);
			float highs = Float.parseFloat(asa[6]);

			if (highs > 4)
				highs = 4;
			if (lows > 4)
				lows = 4;
			Bean b = null;
			// String chr0 = "chr" + asa[1];
			String chr0 = asa[0];
			int pos = Integer.parseInt(asa[2]);

			String key = chr0 + cbhg18.getBand(chr0, pos, pos);
			if (mapArray.containsKey(key)) {
				b = mapArray.get(key);
			} else {
				b = new Bean();
				mapArray.put(key, b);
			}
			// if(chr0.equals("chr1")){
			// System.out.println(as);
			// System.out.println(chr0+"\t"+pos+"\t"+key+"\t"+lows+"\t"+highs);
			// }
			b.high.addValue(highs);
			b.low.addValue(lows);
		}

		return mapArray;

	}

	private static Map<String, Bean> exec2(String id, File arrayf,
			ChromBand cbhg18, Map<String, Float> tp) throws IOException,
			SQLException {

		Map<String, Bean> mapArray = new LinkedHashMap<String, Bean>();

		BufferedReader abr = new BufferedReader(new FileReader(arrayf));
		String as = null;
		abr.readLine();// title

		SummaryStatistics hiave = new SummaryStatistics();
		SummaryStatistics lowave = new SummaryStatistics();

		while ((as = abr.readLine()) != null) {
			String[] asa = as.split("\t");
			// float lows = Float.parseFloat(asa[5]);
			// float highs = Float.parseFloat(asa[6]);

			float lows = Float.parseFloat(asa[6]);
			float highs = Float.parseFloat(asa[5]);

			float tpv = tp.get(id);

			boolean notnoisy = check(id);
			if (notnoisy) {
				lows = adjustlow(lows, tpv);
				highs = adjusthigh(highs, tpv);
			}
			if (highs > 4)
				highs = 4;
			if (lows > 4)
				lows = 4;
			if (!Double.isNaN(highs)) {
				hiave.addValue(highs);
			}
			if (!Double.isNaN(lows)) {
				lowave.addValue(lows);
			}

		}

		abr = new BufferedReader(new FileReader(arrayf));
		abr.readLine();// title

		while ((as = abr.readLine()) != null) {
			String[] asa = as.split("\t");
			// float lows = Float.parseFloat(asa[5]);
			// float highs = Float.parseFloat(asa[6]);

			float lows = Float.parseFloat(asa[6]);
			float highs = Float.parseFloat(asa[5]);

			float tpv = tp.get(id);

			boolean notnoisy = check(id);
			if (notnoisy) {
				lows = adjustlow(lows, tpv);
				highs = adjusthigh(highs, tpv);
			}
			if (highs > 4)
				highs = 4;
			if (lows > 4)
				lows = 4;
			Bean b = null;
			// String chr0 = "chr" + asa[1];
			String chr0 = asa[0];
			int pos = Integer.parseInt(asa[1]);

			String key = chr0 + cbhg18.getBand(chr0, pos, pos);
			if (mapArray.containsKey(key)) {
				b = mapArray.get(key);
			} else {
				b = new Bean();
				mapArray.put(key, b);
			}
			// if(chr0.equals("chr1")){
			// System.out.println(as);
			// System.out.println(chr0+"\t"+pos+"\t"+key+"\t"+lows+"\t"+highs);
			// }
			System.out.println(hiave.getMean());
			System.out.println(lowave.getMean());

			b.high.addValue(highs - (hiave.getMean() - 1));
			b.low.addValue(lows - (lowave.getMean() - 1));

		}

		return mapArray;

	}

	private static boolean check(String id) {

		if (id.contains("18-1"))
			return false;
		if (id.contains("7_2"))
			return false;

		return true;
	}

	private static float adjusthigh(float highs, float tpv) {

		if (highs < 1) {
			return highs;
		}

		float diff = highs - 1;
		float ad = diff / tpv;
		return 1 + ad;
	}

	private static float adjustlow(float lows, float tpv) {

		if (lows > 1) {
			return lows;
		}

		float diff = 1 - lows;
		float ad = diff / tpv;
		return 1 - ad;

	}

	private static File getFile(String array, String id) {

		File f = new File(array);
		for (File ff : f.listFiles()) {
			if (ff.getName().startsWith(id)) {
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
						|| name.contains("cnvAllelicDepth");
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
			if (ff.isFile() && ff.getName().contains("cnvAllelicDepth")) {

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
						|| name.contains("textdata.txt");
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
			if (ff.isFile() && ff.getName().contains("textdata.txt")) {

				// System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

}
