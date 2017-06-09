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

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.bean.BaitSampling;
import jp.ac.utokyo.rcast.karkinos.utils.Interval;

public class Coverage implements Serializable {

	final static int bait_sample_size = 50000;

	class IntervalCov implements Serializable {

		Interval iv;
		int cnt;

		public String getInfoStr() {
			return iv.getChrom() + ":" + iv.getStart() + "-" + iv.getEnd() + "\t" + cnt;
		}

	}

	class Container implements Serializable {
		long cnt = 0;

		void inc(int i) {
			cnt = cnt + i;
		}
	}

	List<IntervalCov> clist = new ArrayList<IntervalCov>();

	public void addBinCount(Interval iv, int cnt) {

		IntervalCov ivc = new IntervalCov();
		ivc.cnt = cnt;
		ivc.iv = iv;
		clist.add(ivc);

	}

	Map<String, Container> readcounts10K = new LinkedHashMap<String, Container>();

	public void setCoverageInfo(SAMRecord sam) {

		String chr = sam.getReferenceName();
		int start = sam.getAlignmentStart();
		int serial = start / 10000;
		String key = chr + "-" + serial;
		if (readcounts10K.containsKey(key)) {
			Container cn = readcounts10K.get(key);
			cn.inc(sam.getReadLength());
		} else {
			Container cn = new Container();
			cn.inc(sam.getReadLength());
			readcounts10K.put(key, cn);
		}

	}

	public Map<String, Container> getReadcounts10K() {
		return readcounts10K;
	}

	public void setReadcounts10K(Map<String, Container> readcounts10k) {
		readcounts10K = readcounts10k;
	}

	Map<CapInterval, Container> capregion = new LinkedHashMap<CapInterval, Container>();
	boolean mapinit = false;
	long maptotal = 0;

	public double millionAdjust() {
		return (double) ((double) 10000000 / (double) maptotal);
	}

	long regcnt = 0;
	CapInterval ciprev = null;

	SAMRecord b4sam = null;
	Map<String, SAMRecord> map = new HashMap<String, SAMRecord>();

	public OntagetInfo setCaptureInfo(SAMRecord sam, CaptureHolder ch) {

		OntagetInfo ret = new OntagetInfo();
		boolean inTarget = false;
		if (mapinit == false) {
			mapinit(ch);
			mapinit = true;
		}

		if (b4sam != null && b4sam.getReferenceIndex() != sam.getReferenceIndex()) {
			// new chrom
			map = null;
			map = new HashMap<String, SAMRecord>();

		}

		CapInterval ci = null;
		if (ciprev != null && ciprev.intersect(sam)) {
			ci = ciprev;
		} else {
			ci = ch.getOverlapping(sam);
		}
		//
		boolean onmarjin = false;
		if (ci == null) {
			ci = ch.getOverlapping(sam, KarkinosProp.baitmergin);
			onmarjin = true;
		}

		if (ci != null) {

			ciprev = ci;
			inTarget = true;
			Container cn = null;
			if (capregion.containsKey(ci)) {
				cn = capregion.get(ci);
			} else {
				cn = new Container();
				capregion.put(ci, cn);
			}
			cn.inc(matchedlength(sam, map));
			regcnt++;
			if (regcnt < bait_sample_size) {
				if (regcnt % 1000 == 0)
					System.out.println(regcnt);
				regBaitSampling(ci, sam);
			}
		}
		if (!inTarget && sam.getReadPairedFlag()) {
			inTarget = (ch.getOverlappingOfPair(sam) != null);
		}
		if (inTarget) {
			maptotal++;
		}
		ret.setOntag(inTarget);
		ret.setOnmarjin(onmarjin);
		b4sam = sam;
		
		return ret;
	}

	private int matchedlength(SAMRecord sam, Map<String, SAMRecord> map) {

		int dlen = mlen(sam); 
		int len = 0;
		if(map.containsKey(sam.getReadName())){
			
			SAMRecord pair = map.get(sam.getReadName());
			//
			if(pair.getAlignmentEnd()>sam.getAlignmentStart()){
				//
				len = sam.getAlignmentEnd()-pair.getAlignmentEnd();
				if(len<0 || len>dlen){
					len = dlen;
				}
				
			}
			map.remove(sam.getReadName());
			
		}else{
			map.put(sam.getReadName(), sam);
		}
		return len;
		
	}

	private int mlen(SAMRecord sam) {

		int ret = sam.getReadLength();
		int ret2 = 0;
		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH)) {

				//
				ret2 = ret2 + ce.getLength();

			}

		}
		if (ret2 > 0)
			return ret2;
		return ret;
	}

	BaitSampling bs = new BaitSampling();

	private void regBaitSampling(CapInterval ci, SAMRecord sam) {

		//
		bs.regBaitSampling(ci, sam);
	}

	public void analysisBaitSampling() {
		bs.analyze();
	}

	public BaitSampling getBs() {
		return bs;
	}

	public void setBs(BaitSampling bs) {
		this.bs = bs;
	}

	private void mapinit(CaptureHolder ch) {

		Iterator<String> ite = ch.map.keySet().iterator();
		while (ite.hasNext()) {

			String key = ite.next();
			Map<Integer, CapInterval> cmap = ch.map.get(key);
			System.out.println(key);
			Iterator<Integer> itepos = cmap.keySet().iterator();
			while (itepos.hasNext()) {
				Integer ikey = itepos.next();
				CapInterval ci = cmap.get(ikey);
				capregion.put(ci, new Container());
			}
		}

	}

	public void merge(Coverage cov2) {

		//
		regcnt = regcnt + cov2.regcnt;
		maptotal = maptotal + cov2.maptotal;

		Set<Entry<CapInterval, Container>> set = cov2.capregion.entrySet();
		for (Entry<CapInterval, Container> e : set) {
			//
			if (capregion.containsKey(e.getKey())) {

				Container c = capregion.get(e.getKey());
				Container c2 = e.getValue();
				c.cnt = c.cnt + c2.cnt;
				capregion.put(e.getKey(), c);

			} else {
				capregion.put(e.getKey(), e.getValue());
			}
		}

	}

}
