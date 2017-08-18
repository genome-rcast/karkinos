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

import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

public class ReadsCounter implements java.io.Serializable {

	int totalmap = 0;
	int totalunique = 0;
	int totalOnTarget = 0;
	int totalUniqueOntarget = 0;
	int duplicatereads = 0;

	int firstReads = 0;
	int secondReads = 0;

	int bothmap = 0;
	int propermap = 0;

	boolean pairStats = false;
	TreeMap<Integer, CounterA> inserSizeMap = new TreeMap<Integer, CounterA>();

	public TreeMap<Integer, CounterA> getInserSizeMap() {
		return inserSizeMap;
	}

	int count4insert = 0;
	double suminsersize = 0;

	public float getMeanInsertSize() {
		if (count4insert == 0)
			return 0f;

		double d = suminsersize / (double) count4insert;
		return (float) d;
	}

	void inserSizeAdd(int insertSize) {

		if (insertSize < -100 || insertSize > 2000) {
			return;
		}
		count4insert++;
		suminsersize = suminsersize + insertSize;
		CounterA cnt = null;
		if (!inserSizeMap.containsKey(insertSize)) {
			cnt = new CounterA();
			inserSizeMap.put(insertSize, cnt);
		} else {
			cnt = inserSizeMap.get(insertSize);
		}
		cnt.inc();

	}

	public int getDuplicatereads() {
		return duplicatereads;
	}

	public int getTotalmap() {
		return totalmap;
	}

	public int getTotalunique() {
		return totalunique;
	}

	public int getTotalOnTarget() {
		return totalOnTarget;
	}

	public int getTotalUniqueOntarget() {
		return totalUniqueOntarget;
	}

	public float onTargetParcent() {
		double r = (double) ((double) totalOnTarget / (double) totalmap);
		r = r * 100;
		return (float) r;
	}

	public void inc(SAMRecord sam, boolean onTarget) {

		if (sam.getDuplicateReadFlag()) {
			duplicatereads++;
			return;
		}

		totalmap++;
		boolean unique = unique(sam);
		if (unique) {
			totalunique++;
		}
		if (onTarget) {
			totalOnTarget++;
			if (unique) {
				totalUniqueOntarget++;
			}

		}

		if (sam.getReadPairedFlag()) {

			pairStats = true;
			if (sam.getProperPairFlag()) {
				propermap++;
			}
			boolean mapped = !sam.getReadUnmappedFlag();
			if (sam.getAlignmentStart() == 0)
				mapped = false;
			boolean matemapped = !sam.getMateUnmappedFlag();
			if (sam.getMateAlignmentStart() == 0)
				matemapped = false;
			boolean bothmapped = mapped && matemapped;
			if (bothmapped) {
				bothmap++;
			}
			// int insertsize =sam.getInferredInsertSize();
			int insertsize = 0;
			if (bothmapped) {

				int start = sam.getAlignmentStart();
				int end = sam.getMateAlignmentStart();
				if (start > end) {
					int tmp = start;
					start = end;
					end = tmp;
				}
				end = end + sam.getReadLength();
				insertsize = end - start;
				inserSizeAdd(insertsize);
			}

			if (sam.getFirstOfPairFlag()) {
				firstReads++;
			} else {
				secondReads++;
			}

		}

	}

	public boolean isPairStats() {
		return pairStats && bothmap > 0;
	}

	public int getFirstReads() {
		return firstReads;
	}

	public int getSecondReads() {
		return secondReads;
	}

	public int getBothmap() {
		return bothmap;
	}

	public int getPropermap() {
		return propermap;
	}

	private boolean unique(SAMRecord sam) {

		Integer nh = null;
		Integer x0 = null;
		try {
			nh = sam.getIntegerAttribute("NH");
			x0 = sam.getIntegerAttribute("X0");
		} catch (Exception ex) {
			return false;
		}

		if (nh == null && x0 == null) {
			// no information to judge
			return true;
		} else {
			if (nh != null) {
				return nh == 1;
			}
			if (x0 != null) {
				return x0 == 1;
			}
		}
		return false;
	}

	public void merge(ReadsCounter rc) {

		totalmap = totalmap + rc.getTotalmap();
		totalunique = totalunique + rc.getTotalunique();
		totalOnTarget = totalOnTarget + rc.getTotalOnTarget();
		totalUniqueOntarget = totalUniqueOntarget + rc.getTotalUniqueOntarget();
		duplicatereads = duplicatereads + rc.getDuplicatereads();

		firstReads = firstReads + rc.getFirstReads();
		secondReads = secondReads + rc.getSecondReads();

		bothmap = bothmap + rc.getBothmap();
		propermap = propermap + rc.getPropermap();
		count4insert = count4insert + rc.count4insert;
		suminsersize = suminsersize + rc.suminsersize;
		//
		try {
			Set<Entry<Integer, CounterA>> set = rc.inserSizeMap.entrySet();
			for (Entry<Integer, CounterA> et : set) {

				int key = et.getKey();
				CounterA value = et.getValue();
				if (inserSizeMap.containsKey(key)) {
					CounterA thisval = inserSizeMap.get(key);
					thisval.cnt = thisval.cnt + value.cnt;
				} else {
					inserSizeMap.put(key, value);
				}

			}
		} catch (Exception ex) {
		}

	}

}
