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
package jp.ac.utokyo.rcast.karkinos.exec;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.utils.Interval;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class CaptureHolder implements java.io.Serializable {

	protected Map<String, TreeMap<Integer, CapInterval>> map = new LinkedHashMap<String, TreeMap<Integer, CapInterval>>();
	long totalcnt = 0;
	long totallen = 0;

	public static void main(String[] arg) {

		CaptureHolder inst = new CaptureHolder();
		//String bed = "/GLUSTER_DIST/data/users/ueda/testframework/script/vcrome2.1.bed";
		String bed = "/GLUSTER_DIST/data/Genomes/karkinos/genome/SureSelectV5plusLincRNA_Regions.bed";
		String tb = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(tb));

		try {
			inst.loadTargetBed(bed, tgr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public SortedMap<Integer, CapInterval> getIntersectCapinterval(Interval iv) {

		//
		TreeMap<Integer, CapInterval> tm = map.get(iv.getChrom());
		int start = iv.getStart();
		int end = iv.getEnd();
		if (tm == null)
			return null;
		SortedMap sm = tm.subMap(start, end);
		return sm;
	}

	public Map<String, TreeMap<Integer, CapInterval>> getMap() {
		return map;
	}

	public void loadTargetBed(String bed, TwoBitGenomeReader tgr)
			throws IOException {

		if (bed == null) {

			System.out.println("no Target File was specified WGS assumed");

			loadUnitTarget(tgr);
			return;
		}

		System.out.println("TargetFile " + bed + " read");
		File f = new File(bed + ".capregion");
		File f0 = new File(bed + ".capwithdepth");

		//
		if (f0.exists()) {
			loadTargetFromCapregionFile(f0, true);
		} else if (f.exists()) {
			loadTargetFromCapregionFile(f, false);
		} else {
			loadTargetBedFirstTime(new File(bed), tgr);
		}
	}

	private void loadUnitTarget(TwoBitGenomeReader tgr) throws IOException {

		// assign 1K artifitial interval
		Map<String, Integer> chrmap = tgr.getReadSizes();
		Iterator<String> ite = chrmap.keySet().iterator();
		int unit  = 10000;
		while (ite.hasNext()) {
			
			String chr = ite.next();
			TreeMap<Integer, CapInterval> tm = map.get(chr);
			if (tm == null) {
				tm = new TreeMap<Integer, CapInterval>();
				map.put(chr, tm);
			}
			//
			int chrend = chrmap.get(chr);
			int start = 1;
			int end = start+unit;
			while(start+unit<chrend){
				
				float cgParcent = tgr.getCGParcent(chr, start, end);
				//debug
				//float cgParcent = 0.5f;
				CapInterval iv = new CapInterval(chr, start, end, false, cgParcent,
						1);			
				tm.put(start, iv);
				start = start + unit;
				end = end + unit;
				
			}		
			

		}

	}

	private void loadTargetFromCapregionFile(File f, boolean withdepth)
			throws IOException {

		BufferedReader br = new BufferedReader(new InputStreamReader(
				new FileInputStream(f)));

		long totallengene = 0;
		long totallenrna = 0;

		try {

			for (;;) {
				String line = br.readLine();
				if (line == null)
					break;
				totalcnt++;
				String[] sa = line.split("\t");
				String chr = sa[0];

				int start = Integer.parseInt(sa[1]);
				int end = Integer.parseInt(sa[2]);
				int length = end - start;
				float cgp = Float.parseFloat(sa[4]);
				float duality = Float.parseFloat(sa[5]);
				boolean gene = sa[6].equals("gene");
				CapInterval iv = new CapInterval(chr, start, end, gene, cgp,
						duality);
				TreeMap<Integer, CapInterval> tm = map.get(chr);
				if (tm == null) {
					tm = new TreeMap<Integer, CapInterval>();
					map.put(chr, tm);
				}
				totallen = totallen + Math.abs(length);
				if (withdepth) {
					float normaldepth = Float.parseFloat(sa[7]);
					iv.setNormalAveDepth(normaldepth);
				}
				if (gene) {
					totallengene = totallengene + Math.abs(length);
				} else {
					totallenrna = totallenrna + Math.abs(length);
				}

				tm.put(start, iv);
			}

		} finally {
			br.close();
		}
		System.out.println("target rigion with " + totalcnt + " region and "
				+ totallen + " bp has loaded");
		System.out.println("gene " + totallengene + "bp has loaded");
		System.out.println("ncRNA " + totallenrna + "bp has loaded");
	}

	public long getTotalcnt() {
		return totalcnt;
	}

	public long getTotallen() {
		return totallen;
	}

	public void loadTargetBedFirstTime(File bed, TwoBitGenomeReader tgr)
			throws IOException {

		// load bed
		BufferedReader br = new BufferedReader(new InputStreamReader(
				new FileInputStream(bed)));
		List<CapInterval> list = new ArrayList<CapInterval>();
		int totalcnt = 0;
		long totallen = 0;
		try {

			for (;;) {
				String line = br.readLine();
				if (line == null)
					break;
				totalcnt++;

				String[] sa = line.split("\t");
				boolean gene = sa.length > 4;

				String chr = sa[0];
				boolean chrContain = tgr.isRefExsist("chr1");
				if (chrContain) {
					if (!chr.contains("chr")) {
						chr = "chr" + chr;
					}
				}
				int start = Integer.parseInt(sa[1]);
				int end = Integer.parseInt(sa[2]);
				CapInterval iv = new CapInterval(chr, start, end, gene);
				list.add(iv);
				totallen = totallen + Math.abs(end - start);
			}

		} finally {
			br.close();
		}
		System.out.println("target rigion with " + totalcnt + " region and "
				+ +totallen + " bp has loaded");
		Collections.sort(list, new MYComparator());
		for (CapInterval iv : list) {

			TreeMap<Integer, CapInterval> tm = map.get(iv.chr);
			if (tm == null) {
				tm = new TreeMap<Integer, CapInterval>();
				map.put(iv.chr, tm);
			}
			Integer before = tm.floorKey(iv.start);
			boolean merge = false;
			if (before != null) {
				CapInterval ivb4 = tm.get(before);
				if (ivb4.intersect(iv)) {
					ivb4.merge(iv);
					merge = true;
				}
			}
			if (merge == false) {
				tm.put(iv.start, iv);
			}

		}

		File f = new File(bed + ".capregion");
		Iterator<String> ite = map.keySet().iterator();
		try {
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(f)));

			while (ite.hasNext()) {

				TreeMap<Integer, CapInterval> tm = map.get(ite.next());
				Set<Entry<Integer, CapInterval>> es = tm.entrySet();
				for (Entry<Integer, CapInterval> entry : es) {

					CapInterval cp = entry.getValue();
					String chr = cp.getChr();
					int start = cp.getStart();
					int end = cp.getEnd();
					int length = end - start;
					float cgParcent = tgr.getCGParcent(chr, start, end);
					float duality = cp.getDuality();
					cp.setCgParcent(cgParcent);
					String genestr = "rna";
					if (cp.gene) {
						genestr = "gene";
					}
					bw.write(chr + "\t" + start + "\t" + end + "\t" + length
							+ "\t" + cgParcent + "\t" + duality + "\t"
							+ genestr + "\n");

				}

			}
			bw.close();
		} catch (Exception ex) {
		}

	}

	class MYComparator implements Comparator<CapInterval> {

		public int compare(CapInterval o1, CapInterval o2) {
			if (o1.chr.equals(o2.chr)) {
				return o1.start - o2.start;
			} else {

				String chr1 = o1.chr;
				if (o1.chr.contains("chr")) {
					chr1 = o1.chr.replace("chr", "");
				}

				String chr2 = o2.chr;
				if (o2.chr.contains("chr")) {
					chr2 = o2.chr.replace("chr", "");
				}
				if (number(chr1) && number(chr2)) {
					return toInt(chr1) - toInt(chr2);
				} else if (number(chr1)) {
					return -1;
				} else if (number(chr2)) {
					return 1;
				} else {
					return chr1.compareTo(chr2);
				}

			}
		}

		private int toInt(String chr) {
			try {
				return Integer.parseInt(chr);
			} catch (NumberFormatException nfe) {
				return -1;
			}
		}

		private boolean number(String chr) {
			try {
				Integer.parseInt(chr);
			} catch (NumberFormatException nfe) {
				return false;
			}
			return true;
		}

	}

	public CapInterval getCapInterval(String chr, int pos) {

		TreeMap<Integer, CapInterval> refmap = map.get(chr);
		if (refmap == null)
			return null;
		Entry<Integer, CapInterval> et = refmap.floorEntry(pos);
		if (et == null) {
			return null;
		} else {

			if (et.getValue().include(pos)) {
				return et.getValue();
			}

		}
		return null;

	}

	public CapInterval getOverlapping(SAMRecord sam) {

		int s = sam.getAlignmentStart();
		int e = sam.getAlignmentEnd();
		if (e == 0 || e == s) {
			e = s + sam.getReadLength();
		}
		return getOverlapping(sam.getReferenceName(), s, e);

	}

	public CapInterval getOverlappingOfPair(SAMRecord sam) {

		if (!sam.getReadPairedFlag()) {
			return null;
		}
		// check only proper
		if (sam.getProperPairFlag()) {
			int s = sam.getMateAlignmentStart();
			int e = s + sam.getReadLength();
			return getOverlapping(sam.getReferenceName(), s, e);
		}
		return null;
	}

	public CapInterval getOverlapping(String chrom, int start, int end) {

		TreeMap<Integer, CapInterval> refmap = map.get(chrom);
		if (refmap == null)
			return null;
		Entry<Integer, CapInterval> et = refmap.floorEntry(end);
		if (et == null) {
			return null;
		} else {

			if (et.getValue().intersect(start, end)) {
				return et.getValue();
			}
		}

		return null;

	}

	public CapInterval getOverlapping(SAMRecord sam, int baitmergin) {
		if (!sam.getReadPairedFlag()) {
			return null;
		}
		// check only proper
		if (sam.getProperPairFlag()) {
			int s = sam.getMateAlignmentStart();
			int e = s + sam.getReadLength();
			return getOverlapping(sam.getReferenceName(), s - baitmergin, e
					+ baitmergin);
		}
		return null;
	}

}
