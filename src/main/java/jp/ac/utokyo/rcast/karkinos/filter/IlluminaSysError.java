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
package jp.ac.utokyo.rcast.karkinos.filter;

import java.io.File;
import java.io.IOException;

import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class IlluminaSysError {

	public static void main(String[] arg) throws IOException {

		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		// String str
		// ="ATAAATTTAAAAGGGAAACTAATTTGGAAATCAGAAAACCACTAAGGAATTTGGGAATTAGGCTTCTGCTGCCCTCTCTGC";
		// boolean b = checkTypicalError('A','C',str);
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		String chrom = "chr1";
		int pos = 155287631;
		String before10 = tgr.getGenomicSeq(chrom, pos - 11, pos - 1, true);
		String after10 = tgr.getGenomicSeq(chrom, pos + 1, pos + 11, true);
		System.out.println(before10 + "\t" + after10);
		String str = tgr.getGenomicSeq(chrom, pos - 40, pos + 40, true);

		String before2 = getBf(before10, 2);
		String after2 = getAfter(after10, 2);
		System.out.println(before2 + "\t" + after2);

		int polynuc = polynuc('G', before10, after10);
		System.out.println(polynuc);
		
		boolean pa1 = polyAorT(before10,6);
		boolean pa2 = polyAorT(after10,6);

		if (pa1 || pa2) {
			System.out.println("true");
		}

		// boolean b = checkTypicalError('A','C',str);

	}

	public static boolean checkTypicalError(char genomeR, char altTumor,
			String neighbor60, int numsupportreads) {

		neighbor60 = neighbor60.toUpperCase();
		String before30 = neighbor60.substring(0, 30);
		String after30 = neighbor60.substring(31);

		String before10 = neighbor60.substring(20, 30);
		String after10 = neighbor60.substring(31, 41);

		// System.out.println(before10+"\t"+after10);

		String before1 = getBf(before30, 1);
		String after1 = getAfter(after30, 1);

		String before2 = getBf(before30, 2);
		String after2 = getAfter(after30, 2);

		String before3 = getBf(before30, 3);
		String after3 = getAfter(before30, 3);

		// String mid20 = neighbor80.substring(30, 50);
		int repeatthres = 12;
		boolean lowsupport = numsupportreads <= 5;
		// boolean lowlowsupport = numsupportreads <= 3;
		int polynuc = polynuc(genomeR, before10, after10);

		if (lowsupport) {
			repeatthres = 10;
			if (polynuc >= 4) {
				return true;
			}
			boolean samenuc = sameNuc(before2 + after1, genomeR, altTumor);
			if (samenuc) {
				return true;
			}
			boolean samenuc2 = sameNuc(before1 + after2, genomeR, altTumor);
			if (samenuc2) {
				return true;
			}

			boolean containnuc = containNuc(before1 +after1,genomeR, altTumor);
			boolean samenucSingle = sameNuc(before1 +after1,genomeR, altTumor);
			
			if (genomeR == 'C' && altTumor == 'A') {
				
				if(containnuc&&numsupportreads==4){
					return true;
				}
				if(samenucSingle){
					return true;
				}
				
			}

			if (genomeR == 'G' && altTumor == 'T') {
				
				if(containnuc&&numsupportreads==4){
					return true;
				}
				if(samenucSingle){
					return true;
				}
				
			}
			
			boolean pa1 = polyAorT(before10,5);
			boolean pa2 = polyAorT(after10,5);
			if (pa1 || pa2) {
				return true;
			}

		}

		if (polynuc >= 6) {
			return true;
		}

		if (genomeR == 'C' && altTumor == 'T') {

			if (before2.equals("GG") || after2.equals("CC")) {
				return true;
			}

		}

		if (genomeR == 'G' && altTumor == 'T') {

			int countGGC = countCodon(before30, "GGC");
			if (before2.equals("GG") || countGGC > 1) {
				return true;
			}
			if (after2.equals("GG")) {
				return true;
			}
			if (before3.equals("TTT")) {
				return true;
			}
			if (after3.equals("GGG")) {
				return true;
			}

			if (moreThanRatio(before10, 'G', 'T', 0.8)) {
				return true;
			}
			if (moreThanRatio(after10, 'G', 'T', 0.8)) {
				return true;
			}
			if (lowsupport) {
				return true;
			}
			repeatthres = 8;

		}

		if (genomeR == 'C' && altTumor == 'A') {

			int countCGG = countCodon(after30, "CGG");
			if (before2.equals("CC")) {
				return true;
			}
			if (after2.equals("CC") || countCGG > 1) {
				return true;
			}
			if (moreThanRatio(before10, 'C', 'A', 0.8)) {
				return true;
			}
			if (moreThanRatio(after10, 'C', 'A', 0.8)) {
				return true;
			}
			if (lowsupport) {
				return true;
			}
			repeatthres = 8;
		}

		boolean pa1 = polyAorT(before10,6);
		boolean pa2 = polyAorT(after10,6);
		if (pa1 || pa2) {
			return true;
		}

		if (contaionsInvertedRepeat(neighbor60, repeatthres)) {
			return true;
		}

		return false;
	}

	private static boolean sameNuc(String s, char genomeR, char altTumor) {

		for (char c : s.toUpperCase().toCharArray()) {

			if (c == genomeR || c == altTumor) {

			} else {
				return false;
			}
		}
		return true;
	}

	private static boolean containNuc(String s, char genomeR, char altTumor) {

		int n= 0;
		for (char c : s.toUpperCase().toCharArray()) {

			if (c == genomeR || c == altTumor) {
				n++;
			} 
		}
		return n>0;
	}

	private static int polynuc(char genomeR, String before, String after) {

		int cnt = 1;
		before = before.toUpperCase();
		after = after.toUpperCase();

		for (int n = before.length() - 1; n >= 0; n--) {
			char c = before.charAt(n);
			if (c == genomeR) {
				cnt++;
			} else {
				break;
			}
		}
		for (char c : after.toCharArray()) {

			//
			if (c == genomeR) {
				cnt++;
			} else {
				break;
			}
		}
		return cnt;
	}

	private static boolean polyAorT(String s,int n) {

		int macp = Math.max(poly(s, 'A'), poly(s, 'T'));
		return macp >= n;

	}

	private static int poly(String s, char ref) {

		int cntmax = 0;
		int cnt = 0;
		for (char c : s.toUpperCase().toCharArray()) {

			//
			if (c == ref) {
				cnt++;
				if (cnt > cntmax) {
					cntmax = cnt;
				}
			} else {
				cnt = 0;
			}

		}
		return cntmax;
	}

	private static boolean moreThanRatio(String mid20, char c, char d, double e) {

		int total = 0;
		int countmatch = 0;
		for (char ca : mid20.toCharArray()) {
			//
			total++;
			if (ca == c || ca == d) {

				//
				countmatch++;
			}

		}
		//
		return ((double) countmatch / (double) total) >= e;

	}

	private static int countCodon(String seq, String s) {
		//
		int cnt = 0;
		for (int n = 0; n < seq.length() - 3; n++) {

			String sseq = seq.substring(n, n + 3);
			if (sseq.equals(s)) {
				cnt++;
			}
		}
		return cnt;
	}

	private static boolean contaionsInvertedRepeat(String ref, int thres) {

		//
		String revcon = CalcUtils.revcon(ref);

		StringWithIndex refwi = new StringWithIndex(ref);
		StringWithIndex revconwi = new StringWithIndex(revcon);

		for (int n = 0; n < ref.length() - thres; n++) {

			refwi.setStartEnd(0, ref.length() - n);
			revconwi.setStartEnd(n, ref.length());
			int maxmatch = maxmatch(refwi, revconwi);
			if (maxmatch >= thres) {
				return true;
			}

		}
		return false;

	}

	private static int maxmatch(StringWithIndex refwi, StringWithIndex revconwi) {

		int loop = refwi.getEnd() - refwi.getStart();
		int maxmatch = 0;
		int count = 0;
		for (int n = 0; n < loop; n++) {

			//
			char refc = refwi.getCharAt(n);
			char revc = revconwi.getCharAt(n);
			// System.out.println(refc+"-"+revc);
			//
			if (refc == revc) {

				//
				count++;
				if (count > maxmatch) {
					maxmatch = count;
				}
			} else {
				count = 0;
			}

		}
		return maxmatch;
	}

	private static String getAfter(String s, int n) {

		return s.substring(0, n);

	}

	private static String getBf(String s, int n) {
		if (s.length() >= n)
			return s.substring(s.length() - n);
		return "";
	}

}
