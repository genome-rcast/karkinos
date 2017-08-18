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
package jp.ac.utokyo.rcast.karkinos.bean;

import htsjdk.samtools.SAMRecord;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class BaitSampling implements Serializable {

	public SummaryStatistics[] getS_fstat() {
		return s_fstat;
	}

	public SummaryStatistics[] getS_rstat() {
		return s_rstat;
	}

	public SummaryStatistics[] getE_fstat() {
		return e_fstat;
	}

	public SummaryStatistics[] getE_rstat() {
		return e_rstat;
	}

	Map<CapInterval, BailSamplingBean> map = new HashMap<CapInterval, BailSamplingBean>();
	public final static int bait_sample_length = 200;
	SummaryStatistics[] s_fstat = new SummaryStatistics[bait_sample_length];
	SummaryStatistics[] s_rstat = new SummaryStatistics[bait_sample_length];

	SummaryStatistics[] e_fstat = new SummaryStatistics[bait_sample_length];
	SummaryStatistics[] e_rstat = new SummaryStatistics[bait_sample_length];

	private void initSS() {

		for (int n = 0; n < bait_sample_length; n++) {
			//
			s_fstat[n] = new SummaryStatistics();
			s_rstat[n] = new SummaryStatistics();
			e_fstat[n] = new SummaryStatistics();
			e_rstat[n] = new SummaryStatistics();

		}

	}

	public void analyze() {

		Set<Entry<CapInterval, BailSamplingBean>> set = map.entrySet();
		initSS();
		for (Entry<CapInterval, BailSamplingBean> e : set) {
			
			int cstart = e.getKey().getStart();
			int cend = e.getKey().getEnd();
			int len = cend-cstart;
			
			System.out.println(cstart+"\t"+cend+"\t"+(cend-cstart));
			BailSamplingBean sample = e.getValue();
			for (int n = 0; n < bait_sample_length; n++) {

				int f1 = sample.fromStartForward[n];
				int r1 = sample.fromStartReverse[n];

				int total = f1 + r1;
				double frate = (double) f1 / (double) total;
				double rrate = (double) r1 / (double) total;
				//
				if (total > 0) {
					s_fstat[n].addValue(frate);
					s_rstat[n].addValue(rrate);
				}

				int f2 = sample.fromEndForward[n];
				int r2 = sample.fromEndReverse[n];
				int total2 = f2 + r2;
				double frate2 = (double) f2 / (double) total2;
				double rrate2 = (double) r2 / (double) total2;
				if (total2 > 0) {
					e_fstat[n].addValue(frate2);
					e_rstat[n].addValue(rrate2);
				}
			}

		}
		//System.out.println("here");

	}

	public void regBaitSampling(CapInterval ci, SAMRecord sam) {

		//
		BailSamplingBean bsb = null;
		if (map.containsKey(ci)) {
			bsb = map.get(ci);
		} else {
			bsb = new BailSamplingBean();
			map.put(ci, bsb);
		}

		// //
		boolean mapNegative = sam.getReadNegativeStrandFlag();

		int rstart = sam.getAlignmentStart();
		int rend = sam.getAlignmentEnd();
		if (rend == 0) {
			rend = rstart + sam.getReadLength();
		}
		// /
		int cstart = ci.getStart();
		int cend = ci.getEnd();

		//
		int relstart = rstart - cstart;
		int relend = rend - cstart;
		if (relstart < 0)
			relstart = 0;
		//
		int erelstart = cend - rend;
		if (erelstart < 0)
			erelstart = 0;
		int erelend = cend - rstart;

		//
		bsb.setForward(relstart, relend, mapNegative);
		bsb.setReverse(erelstart, erelend, mapNegative);

	}

}
