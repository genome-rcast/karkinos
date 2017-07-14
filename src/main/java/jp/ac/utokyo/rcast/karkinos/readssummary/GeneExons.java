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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.utils.DataHolder;
import au.com.bytecode.opencsv.CSVReader;

public class GeneExons implements java.io.Serializable {
	public GeneExons(String refflat) {

		//
		if (refflat != null) {
			loadmap(refflat);
		}

	}
	
	Map<String,DataHolder> counterForGene = new LinkedHashMap<String,DataHolder>();
	long normaltotal;
	long tumortotal;
	

	public Map<String, DataHolder> getCounterForGene() {
		return counterForGene;
	}


	public void setCounterForGene(Map<String, DataHolder> counterForGene) {
		this.counterForGene = counterForGene;
	}


	public long getNormaltotal() {
		return normaltotal;
	}


	public void setNormaltotal(long normaltotal2) {
		this.normaltotal = normaltotal2;
	}


	public long getTumortotal() {
		return tumortotal;
	}


	public void setTumortotal(long tumortotal) {
		this.tumortotal = tumortotal;
	}

	Map<String,Integer> counter = new HashMap<String,Integer>();
	private void loadmap(String refflat) {
		// TODO Auto-generated method stub
		try {
			CSVReader brcvs = new CSVReader(new FileReader(refflat), '\t');

			String[] data = null;
			while ((data = brcvs.readNext()) != null) {

				try {
					
					String geneSymbol = data[0];
					String symbolWithoutNum = symbolWithoutNum(geneSymbol);
					//System.out.println(geneSymbol+"\t"+symbolWithoutNum);
					if(counter.containsKey(symbolWithoutNum)){
						int countn = counter.get(symbolWithoutNum);
						countn++;
						counter.put(symbolWithoutNum, countn);
					}else{
						counter.put(symbolWithoutNum, 1);
					}
					
					
					String refseqid = data[1];
					String chrom = data[2];
					int txstart = Integer.parseInt(data[4]);
					int txend = Integer.parseInt(data[5]);
					TreeMap<Integer, Interval> tm0 = genemap.get(chrom);
					if (tm0 == null) {
						tm0 = new TreeMap<Integer, Interval>();
						genemap.put(chrom, tm0);
					}
					Interval iv0 = new Interval(chrom, txstart, txend, refseqid,geneSymbol);
					tm0.put(txstart, iv0);

					int cdsstart = Integer.parseInt(data[6]);
					int cdsend = Integer.parseInt(data[7]);
					String exonstarts = data[9];
					String exonends = data[10];
					List<Integer> es = toIntList(exonstarts);
					List<Integer> ee = toIntList(exonends);
					// //
					int idx = 0;
					for (int start : es) {
						int end = ee.get(idx);
						if (start < cdsstart && end < cdsstart) {
							continue;
						} else if (start < cdsstart && end >= cdsstart) {
							start = cdsstart;
						} else if (start > cdsend && end > cdsend) {
							continue;
						} else if (start <= cdsend && end > cdsend) {
							end = cdsend;
						}
						// /
						TreeMap<Integer, Interval> tm = map.get(chrom);
						if (tm == null) {
							tm = new TreeMap<Integer, Interval>();
							map.put(chrom, tm);
						}
						//
						Interval iv = new Interval(chrom, start, end);
						iv.end = end;
						iv.geneSymbol = geneSymbol;
						iv.refseqid = refseqid;
						iv.exonidx = idx+1;
						
						 
						tm.put(start, iv);
						idx++;
					}

				} catch (Exception e) {

				}

			}
			
			

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public String symbolWithoutNum(String s){
		
		StringBuffer sb  = new StringBuffer();
		for(char c :s.toCharArray()){
			
			if(!isNumber(c)){
			 sb.append(c);
			}
		}
		return sb.toString();
		
	}

	private boolean isNumber(char c) {
		if(c<'0'){
			return false;
		}
		if(c>'9'){
			return false;
		}
		return true;
	}


	private List<Integer> toIntList(String s) {

		ArrayList<Integer> list = new ArrayList<Integer>();
		for (String ss : s.split(",")) {
			try {
				int i = Integer.parseInt(ss);
				list.add(i);
			} catch (Exception ex) {
			}
		}
		return list;
	}

	Map<String, TreeMap<Integer, Interval>> map = new LinkedHashMap<String, TreeMap<Integer, Interval>>();

	Map<String, TreeMap<Integer, Interval>> genemap = new LinkedHashMap<String, TreeMap<Integer, Interval>>();

	Interval ivprev;
	
	Interval ivprevExon;
	

	public boolean onCDS(String chr, int pos) {

		if (ivprev != null && ivprev.contain(chr, pos)) {
			return true;
		} else {
			Interval iv = getIV(chr, pos,map);
			if (iv != null) {

				ivprev = iv;
				return iv.contain(chr, pos);

			}
		}

		return false;
	}

	Interval ivprevgene;
	public String getGeneId(String chr, int pos) {

		if (ivprevgene != null && ivprevgene.contain(chr, pos)) {
			return ivprevgene.refseqid;
		} else {
			Interval iv = getIV(chr, pos,genemap);
			if (iv != null) {

				ivprevgene = iv;
				if(iv.contain(chr, pos)){
					return iv.refseqid;
				}

			}
		}

		return null;
	}
	
	public Interval getGeneInterval(String chr, int pos) {

		if (ivprevgene != null && ivprevgene.contain(chr, pos)) {
			return ivprevgene;
		} else {
			Interval iv = getIV(chr, pos,genemap);
			if (iv != null) {

				ivprevgene = iv;
				if(iv.contain(chr, pos)){
					return iv;
				}

			}
		}

		return null;
	}
	
	
	public Interval getGeneIntervalExon(String chr, int pos) {

		if (ivprevExon != null && ivprevExon.contain(chr, pos)) {
			return ivprevExon;
		} else {
			Interval iv = getIV(chr, pos,map);
			if (iv != null) {

				ivprevExon = iv;
				if(iv.contain(chr, pos)){
					return iv;
				}

			}
		}

		return null;
	}
	
	
	public boolean isFrequentIsoform(String chr, int pos){
		
		String gs = getGeneSymbol(chr,pos);
		if(gs==null)return false;
		String gswon = symbolWithoutNum(gs);
		if(counter.containsKey(gswon)){
			int count = counter.get(gswon);
			return count > 200;
		}
		return false;
		
	}
	
	public String getGeneSymbol(String chr, int pos) {

		if (ivprevgene != null && ivprevgene.contain(chr, pos)) {
			return ivprevgene.geneSymbol;
		} else {
			Interval iv = getIV(chr, pos,genemap);
			if (iv != null) {

				ivprevgene = iv;
				if(iv.contain(chr, pos)){
					return iv.geneSymbol;
				}

			}
		}

		return null;
	}
	
	
	

	private Interval getIV(String chr, int pos,Map<String, TreeMap<Integer, Interval>> map) {

		TreeMap<Integer, Interval> tm = map.get(chr);
		if (tm != null) {

			Entry<Integer, Interval> et = tm.floorEntry(pos);
			if (et != null) {
				Interval iv = et.getValue();
				if (iv != null && iv.contain(chr, pos)) {
					return iv;
				}
			}

		}
		return null;
	}


	public String getGeneSymbols(String chr, int start, int end) {
		
		Set<String> s = new LinkedHashSet<String>();
		
		int n = start -1000;
		while(n<end+1000){
			String gs = getGeneSymbol(chr,n);
			if(gs!=null){
				s.add(gs);
			}
			n=n+1000;
		}
		StringBuffer sb = new StringBuffer();
		for(String ss:s){
			sb.append(ss+",");
		}
		return sb.toString();
	}

}
