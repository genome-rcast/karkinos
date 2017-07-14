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

import java.util.LinkedHashMap;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.readssummary.SNPDepthCounter;

public class NormalSNPCounter {

//	String keys[] = new String[] { "AtoT", "AtoG", "AtoC", "TtoA", "TtoG",
//			"TtoC", "CtoA", "CtoT", "CtoG", "GtoA", "GtoT", "GtoC" };

	Map<String, SNPDepthCounter> mutationCounterNormalSNP = new LinkedHashMap<String, SNPDepthCounter>();
	int total = 0;
	public NormalSNPCounter(DataSet dataset) {
		for (String key : GenotypeKeyUtils.keys1) {
			mutationCounterNormalSNP.put(key, new SNPDepthCounter());
		}
		setData(dataset);
	}

	public void setData(DataSet dataset) {

		for (SNVHolder snvSNP : dataset.getSnvlist()) {
			int flg = snvSNP.getFlg();
			if (flg == PileUP.NormalSNP) {

				//
				char ref = snvSNP.getNormal().getGenomeR();
				char alt = snvSNP.getNormal().getALT();
				String key = ref + "to" + alt;
				key = GenotypeKeyUtils.aggrigateKeys(key);
				SNPDepthCounter counter = null;
				boolean hetroSNP = snvSNP.getNormal().getRatio() < 0.6;
				if (hetroSNP && mutationCounterNormalSNP.containsKey(key)) {
					counter = mutationCounterNormalSNP.get(key);
					counter.reg(snvSNP.getNormal().getTotalcnt());
					total++;
				}
			}
		}

	}

	
	public float getHetroSNPRatio(String key){
		//6 allele
		
		SNPDepthCounter counter = null;
		key = GenotypeKeyUtils.aggrigateKeys(key);
		if (mutationCounterNormalSNP.containsKey(key)) {
			
			counter = mutationCounterNormalSNP.get(key);
			
		}
		int cnt = 1;
		if(counter!=null){
			cnt = cnt+counter.getTotal();
		}
		return (float)((double)cnt/(double)(total+6));
		
	}
	public float getHetroSNPRatioRemain(String key){
		return 1-getHetroSNPRatio(key);
		
	}
	
	public float getLogRefHetroSNPRatio(String key, long ntotal){
		
		key = GenotypeKeyUtils.aggrigateKeys(key);
		SNPDepthCounter counter = null;
		if (mutationCounterNormalSNP.containsKey(key)) {
			
			counter = mutationCounterNormalSNP.get(key);
			double d = (double)((double)ntotal/(double)(counter.getTotal()+1));
			return (float)(Math.log10(d));
		}
		return 1;
		
	}
	
}
