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
package jp.ac.utokyo.rcast.karkinos.readssummary;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.stream.Collectors;

import jp.ac.utokyo.rcast.karkinos.utils.DataHolder;
import au.com.bytecode.opencsv.CSVReader;

public class GeneExons implements java.io.Serializable {
	private Map<String,DataHolder> counterForGene = new LinkedHashMap<String,DataHolder>();
	private long normaltotal = 0L;
	private long tumortotal = 0L;
	private Interval ivprev = null;
	private Interval ivprevExon = null; // Cache `ivprev` for exon.
	private Interval ivprevgene = null; // Cache `ivprev`.

	Map<String, Integer> counter = new HashMap<String, Integer>();
	Map<String, TreeMap<Integer, Interval>> map = new LinkedHashMap<String, TreeMap<Integer, Interval>>();
	Map<String, TreeMap<Integer, Interval>> genemap = new LinkedHashMap<String, TreeMap<Integer, Interval>>();

	public GeneExons(final String refflat) throws IOException {
		if (refflat != null) {
			loadmap(refflat);
		}
	}

	public Map<String, DataHolder> getCounterForGene() {
		return counterForGene;
	}

	public long getNormaltotal() {
		return normaltotal;
	}

	public void setNormaltotal(final long normaltotal) {
		this.normaltotal = normaltotal;
	}

	public long getTumortotal() {
		return tumortotal;
	}

	public void setTumortotal(final long tumortotal) {
		this.tumortotal = tumortotal;
	}

	private void loadmap(final String refflat) throws IOException {
		try (final CSVReader brcvs = new CSVReader(new FileReader(refflat), '\t')) {
			String[] data = null;
			while ((data = brcvs.readNext()) != null) {
				final String geneSymbol;
				if (data.length == 11) { // refFlat.txt
					geneSymbol = data[0];
				} else if (data.length == 16) { // refGene.txt
					geneSymbol = data[12];
				} else {
					throw new IOException("Illegal `refFlat` or `refGene` format.");
				}

				final String symbolWithoutNum = symbolWithoutNum(geneSymbol);
				this.counter.put(symbolWithoutNum, 1 + this.counter.getOrDefault(symbolWithoutNum, 0));

				// Convert refGene.txt format (0-based half-opened range)
				// into Interval format (1-based closed open range).
				final String refseqid = data[1];
				final String chrom = data[2];
				final int txStart = Integer.parseInt(data[4]) + 1;
				final int txEnd = Integer.parseInt(data[5]);
				this.genemap
				    .computeIfAbsent(chrom, k -> new TreeMap<>())
				    .put(txStart, new Interval(chrom, txStart, txEnd, refseqid, geneSymbol));

				final int cdsStart = Integer.parseInt(data[6]) + 1;
				final int cdsEnd = Integer.parseInt(data[7]);
				final int exonCount = Integer.parseInt(data[8]);
				final List<Integer> exonStarts = toIntList(data[9]).stream().map(e -> e + 1).collect(Collectors.toList());
				final List<Integer> exonEnds = toIntList(data[10]);

				if (cdsStart <= cdsEnd) {
					for (int i = 0; i < exonCount; ++i) {
						int start = exonStarts.get(i);
						int end = exonEnds.get(i);
						if (end < cdsStart || cdsEnd < start) {
							continue;
						}
						start = Math.max(start, cdsStart);
						end = Math.min(end, cdsEnd);
						final Interval iv = new Interval(chrom, start, end, refseqid, geneSymbol);
						iv.exonidx = i + 1;
						this.map
						    .computeIfAbsent(chrom, k -> new TreeMap<>())
						    .put(start, iv);
					}
				}
			}
		}
	}

	private static String symbolWithoutNum(final String symbol) {
		return symbol.codePoints()
		             .filter(c -> !isNumber((char)c))
		             .mapToObj(c -> String.valueOf((char)c))
		             .collect(Collectors.joining());
	}

	private static boolean isNumber(final char c) {
		return '0' <= c && c <= '9';
	}

	private static List<Integer> toIntList(final String str) {
		final ArrayList<Integer> list = new ArrayList<Integer>();
		for (final String s : str.split(",")) {
			try {
				list.add(Integer.parseInt(s));
			} catch (Exception ex) {
				// FALLTHROUGH
			}
		}
		return list;
	}

	public boolean onCDS(final String chr, final int pos) {
		if (this.ivprev != null && this.ivprev.contain(chr, pos)) {
			return true;
		}
		final Interval iv = this.getIV(chr, pos, this.map);
		if (iv == null) {
			return false;
		}
		this.ivprev = iv;
		return iv.contain(chr, pos);
	}

	public String getGeneId(final String chr, final int pos) {
		if (this.ivprevgene != null && this.ivprevgene.contain(chr, pos)) {
			return this.ivprevgene.refseqid;
		}
		final Interval iv = this.getIV(chr, pos, this.genemap);
		if (iv == null) {
			return null;
		}
		this.ivprevgene = iv;
		return iv.contain(chr, pos) ? iv.refseqid : null;
	}

	public Interval getGeneIntervalExon(final String chr, final int pos) {
		if (this.ivprevExon != null && this.ivprevExon.contain(chr, pos)) {
			return this.ivprevExon;
		}
		final Interval iv = this.getIV(chr, pos, this.map);
		if (iv == null) {
			return null;
		}
		this.ivprevExon = iv;
		return iv.contain(chr, pos) ? iv : null;
	}

	public boolean isFrequentIsoform(final String chr, final int pos) {
		final String s = this.getGeneSymbol(chr, pos);
		if (s == null) return false;
		final String t = this.symbolWithoutNum(s);
		return this.counter.containsKey(t) ? this.counter.get(t) > 200 : false;
	}

	public String getGeneSymbol(final String chr, final int pos) {
		if (this.ivprevgene != null && this.ivprevgene.contain(chr, pos)) {
			return this.ivprevgene.geneSymbol;
		}
		final Interval iv = this.getIV(chr, pos, this.genemap);
		if (iv == null) {
			return null;
		}
		this.ivprevgene = iv;
		return iv.contain(chr, pos) ? iv.geneSymbol : null;
	}

	private static Interval getIV(final String chr, final int pos, final Map<String, TreeMap<Integer, Interval>> map) {
		if (map.containsKey(chr)) {
			final Entry<Integer, Interval> et = map.get(chr).floorEntry(pos);
			if (et != null) {
				final Interval iv = et.getValue();
				if (iv != null && iv.contain(chr, pos)) {
					return iv;
				}
			}
		}
		return null;
	}
}
