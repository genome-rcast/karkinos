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

import htsjdk.samtools.SAMRecord;

import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.utils.Interval;

public class ReadsSummary implements java.io.Serializable {

	public DepthCounter getNormalDepth() {
		return normalDepth;
	}

	public DepthCounter getTumorDepth() {
		return tumorDepth;
	}

	String normalbam;
	String tumorbam;
	String taretbed;

	public String getTaretbed() {
		return taretbed;
	}

	public void setTaretbed(String taretbed) {
		this.taretbed = taretbed;
	}

	public String getNormalbam() {
		return normalbam;
	}

	public void setNormalbam(String normalbam) {
		this.normalbam = normalbam;
	}

	public String getTumorbam() {
		return tumorbam;
	}

	public void setTumorbam(String tumorbam) {
		this.tumorbam = tumorbam;
	}

	ReadsCounter normalCounter = new ReadsCounter();
	ReadsCounter tumorCounter = new ReadsCounter();

	DepthCounter normalDepth = new DepthCounter();
	DepthCounter tumorDepth = new DepthCounter();

	public void setNormalDepth(String chr, int pos, int n, int nontarget) {
		normalDepth.add(chr, pos, n,nontarget);
	}

	public void setTumorDepth(String chr, int pos, int n, int nontarget) {
		tumorDepth.add(chr, pos, n,nontarget);
	}

	public Map<String, ReadsCounter> getNormalPerChrom() {
		return normalPerChrom;
	}

	public Map<String, ReadsCounter> getTumorPerChrom() {
		return tumorPerChrom;
	}

	public ReadsCounter getNormalCounter() {
		return normalCounter;
	}

	public ReadsCounter getTumorCounter() {
		return tumorCounter;
	}

	Set<String> alreadyregset = new HashSet<String>();

	public void resetAlreadyregset() {
		alreadyregset = null;
		alreadyregset = new HashSet<String>();
	}

	Map<String, ReadsCounter> normalPerChrom = new LinkedHashMap<String, ReadsCounter>();
	Map<String, ReadsCounter> tumorPerChrom = new LinkedHashMap<String, ReadsCounter>();

	public void regN(SAMRecord sam, boolean onTarget, Interval iv) {
		//
		if (areadyCounts(sam, iv)) {
			return;
		}

		normalCounter.inc(sam, onTarget);
		String chr = sam.getReferenceName();
		ReadsCounter rc = null;
		if (normalPerChrom.containsKey(chr)) {
			rc = normalPerChrom.get(chr);
		} else {
			rc = new ReadsCounter();
			normalPerChrom.put(chr, rc);
		}
		rc.inc(sam, onTarget);
	}

	private boolean areadyCounts(SAMRecord sam, Interval iv) {

		if (iv == null)
			return false;

		if (iv.getEnd() <= sam.getAlignmentEnd()) {
			alreadyregset.add(key(sam));
		} else if (sam.getAlignmentStart() <= iv.getStart()) {
			return alreadyregset.contains(key(sam));

		}
		return false;
	}

	private String key(SAMRecord sam) {
		String name = sam.getReadName();
		if (sam.getReadPairedFlag()) {
			name = name + sam.getFirstOfPairFlag();
		}
		return name;
	}

	public void regT(SAMRecord sam, boolean onTarget, Interval iv) {
		//
		if (areadyCounts(sam, iv)) {
			return;
		}
		tumorCounter.inc(sam, onTarget);
		String chr = sam.getReferenceName();
		ReadsCounter rc = null;
		if (tumorPerChrom.containsKey(chr)) {
			rc = tumorPerChrom.get(chr);
		} else {
			rc = new ReadsCounter();
			tumorPerChrom.put(chr, rc);
		}
		rc.inc(sam, onTarget);

	}

	public boolean isPairStats() {

		return normalCounter.isPairStats() && tumorCounter.isPairStats();
	}

	public int[] getInsertSizeInterval() {
		int start = 0;
		int startN = 0;
		if (normalCounter.isPairStats()) {
			startN = normalCounter.inserSizeMap.firstKey();
		}
		int startT = 0;
		if (tumorCounter.isPairStats()) {
			startT = tumorCounter.inserSizeMap.firstKey();
		}
		if (startN < 0 || startT < 0) {
			start = Math.min(startN, startT);
		}
		int end = 800;
		int endN = 800;
		if (normalCounter.isPairStats()) {
			endN = normalCounter.inserSizeMap.lastKey();
		}
		int endT = 800;
		if (tumorCounter.isPairStats()) {
			endT = tumorCounter.inserSizeMap.lastKey();
		}
		if (endN > end || endT > end) {
			end = Math.max(endN, endT);
		}
		int[] ary = { start, end };
		return ary;
	}

	public void merge(ReadsSummary rs2) {

		try{
			normalCounter.merge(rs2.getNormalCounter());
			tumorCounter.merge(rs2.getTumorCounter());
		}catch(Exception ex){}
		
		try{
			normalDepth.merge(rs2.getNormalDepth());
			tumorDepth.merge(rs2.getTumorDepth());
		}catch(Exception ex){}
		//
		Iterator<String> ite = rs2.getNormalPerChrom().keySet().iterator();
		while(ite.hasNext()){
			
			String chrom = ite.next();
			ReadsCounter rcn = rs2.getNormalPerChrom().get(chrom);
			if(normalPerChrom.containsKey(chrom)){
				
				ReadsCounter rcno = normalPerChrom.get(chrom);			
				//
				rcno.merge(rcn);				
				
			}else{
				normalPerChrom.put(chrom, rcn);
			
			}
		}		
		
		ite = rs2.getTumorPerChrom().keySet().iterator();
		while(ite.hasNext()){
			
			String chrom = ite.next();
			ReadsCounter rct = rs2.getTumorPerChrom().get(chrom);
			if(tumorPerChrom.containsKey(chrom)){

				ReadsCounter rcto = tumorPerChrom.get(chrom);
				rcto.merge(rct);
				
			}else{
				tumorPerChrom.put(chrom,rct);
			}
		}		
		


	}

	long[] nucrefcount = new long[4];

	public void setNucCountRef(char genomeR) {
		int idx = getNucIndex(genomeR);
		if(idx>=0){
			nucrefcount[idx] = nucrefcount[idx]+1;
		}
	}
	
	public long getNucCountRef(char genomeR) {
		int idx = getNucIndex(genomeR);
		if(idx>=0){
			return nucrefcount[idx];
		}
		return 0;
	}

	long[] nuccount = new long[4];
	public void setNucCount(char genomeR) {
		int idx = getNucIndex(genomeR);
		if(idx>=0){
			nuccount[idx] = nuccount[idx]+1;
		}
	}

	public int getNucIndex(char genomeR) {
		switch (genomeR) {

		case 'A':
			return 0;
		case 'a':
			return 0;
		case 'T':
			return 1;
		case 't':
			return 1;
		case 'G':
			return 2;
		case 'g':
			return 2;
		case 'C':
			return 3;
		case 'c':
			return 3;

		}
		return -1;
	}
	int readslenn = 0;
	public void setReadslenn(int readslen) {
		if(readslen>this.readslenn){
			this.readslenn = readslen;
		}
	}
	int readslent = 0;
	public void setReadslent(int readslen) {
		if(readslen>this.readslent){
			this.readslent = readslen;
		}
	}

	public int getReadslen() {
		return Math.max(readslent,readslenn);
	}
	
	String refflat = null;
	public void setRefFlat(String refflat) {
		this.refflat = refflat;	
		GeneExons ge = new GeneExons(refflat);
		normalDepth.setGeneExons(ge);
		tumorDepth.setGeneExons(ge);
		
	}
	public void clearGeneExon() {
		this.refflat = null;	
		normalDepth.setGeneExons(null);
		tumorDepth.setGeneExons(null);
		
	}
	//
	
	
	

}
