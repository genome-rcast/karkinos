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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CaptureHolder;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;

public class MakeAveDepth extends ReadWriteBase {

	public static void main(String[] arg){
		
		CaptureHolder inst = new CaptureHolder();
		String bed = "/GLUSTER_DIST/data/Genomes/karkinos/genome/halo_ver2_Regions.bed";
		String tb = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(tb));
		try {
			inst.loadTargetBed(bed, tgr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void _main(String[] args) {

		//
		CaptureHolder inst = new CaptureHolder();
		String bed = "/GLUSTER_DIST/data/users/ueda/testframework/script/vcrome2.1.bed";
		String tb = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		String bams = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_Baylor";
		String regexp = ".*normal_genome.bam";
		int mergin = KarkinosProp.baitmergin;
		// String bed = args[0];
		// String tb= args[1];
		// //
		// String bams = args[2];
		// String regexp = args[3];

		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(tb));
		try {
			inst.loadTargetBed(bed, tgr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//
		List<String> list = new ArrayList<String>();
		searchRecursive(new File(bams), list, regexp);
		Map<CapInterval, Double> map = new LinkedHashMap<CapInterval, Double>();
		double total = 0;
		int readlen = 100;
		int samplecnt = 0;
		for (String s : list) {

			if (samplecnt >= 30)
				break;
			if (s.contains("135") || s.contains("183")|| s.contains("186")|| s.contains("324")|| s.contains("412")) {
				continue;
			}

			samplecnt++;

			SAMFileReader normalbamr = getReader(s);

			Set<String> ite = inst.getMap().keySet();
			for (String chr : ite) {

				System.out.println(s + "\t" + chr);
				CloseableIterator<SAMRecord> normarIte = normalbamr.query(chr,
						0, 0, true);

				CapInterval ci = null;
				while (normarIte.hasNext()) {
					
					SAMRecord sam = null;
					try {
						sam = normarIte.next();
					} catch (Exception ex) {
						continue;
					}
					if (sam.getReadLength() != readlen) {
						readlen = sam.getReadLength();
					}
					if (ci == null
							|| !ci.intersect(sam.getAlignmentStart() - mergin,
									sam.getAlignmentEnd() + mergin)) {

						ci = inst.getOverlapping(sam.getReferenceName(),
								sam.getAlignmentStart() - mergin,
								sam.getAlignmentEnd() + mergin);
					}
					if (ci != null) {

						total = total + sam.getReadLength();
						double d = 0;
						if (map.containsKey(ci)) {

							//
							d = map.get(ci);

						}
						d = d + sam.getReadLength();
						map.put(ci, d);

					}

				}
				normarIte.close();
			}
			normalbamr.close();
			System.gc();
		}
		System.out.println(total);

		File f = new File(bed + ".capwithdepth");
		Iterator<String> ite = inst.getMap().keySet().iterator();
		try {
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(f)));

			while (ite.hasNext()) {

				TreeMap<Integer, CapInterval> tm = inst.getMap()
						.get(ite.next());
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
					if (cp.isGene()) {
						genestr = "gene";
					}
					float millionavedepth = 0;
					if (map.containsKey(cp)) {

						double totalbase = map.get(cp);
						double r = totalbase / total;
						r = r * 1000000;
						r = r / cp.getLength();
						millionavedepth = (float) r;
					}
					bw.write(chr + "\t" + start + "\t" + end + "\t" + length
							+ "\t" + cgParcent + "\t" + duality + "\t"
							+ genestr + "\t" + millionavedepth + "\n");

				}

			}
			bw.close();
		} catch (Exception ex) {
		}

	}

	private static boolean searchRecursive(File f, final List<String> list,
			final String reg) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory() || name.matches(reg);

				if (!target)
					return false;
				//
				String absPath = dir.getAbsolutePath() + File.separator + name;
				if (new File(absPath).isFile()) {
					return true;
				} else {
					return searchRecursive(new File(absPath), list, reg);
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

			if (ff.isDirectory())
				continue;
			if (ff.getName().matches(reg)) {

				list.add(ff.getAbsolutePath());

			}
		}
		return true;

	}

}
