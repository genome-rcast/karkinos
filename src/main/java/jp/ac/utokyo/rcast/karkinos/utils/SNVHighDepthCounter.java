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
import jp.ac.utokyo.rcast.karkinos.filter.FilterResult;
import jp.ac.utokyo.rcast.karkinos.readssummary.SNPDepthCounter;

public class SNVHighDepthCounter {

//	String keys[] = new String[] { "AtoT", "AtoG", "AtoC", "TtoA", "TtoG",
//			"TtoC", "CtoA", "CtoT", "CtoG", "GtoA", "GtoT", "GtoC" };

	Map<String, SNPDepthCounter> mutationCounter = new LinkedHashMap<String, SNPDepthCounter>();
	int total = 0;
	public SNVHighDepthCounter(DataSet dataset) {
		for (String key : GenotypeKeyUtils.keys1) {
			
			SNPDepthCounter counter = new SNPDepthCounter();
			counter.reg(100);
			mutationCounter.put(key,counter);
		}
		setData(dataset);
	}

	public void setData(DataSet dataset) {

		for (SNVHolder snv : dataset.getSnvlist()) {
			
			//
			FilterResult fr = snv.getFilterResult();
			if(fr!=null && fr.isPassFilter()){
				
				if(snv.getTumor().getTotalcnt()>100){
					
					char ref = snv.getNormal().getGenomeR();
					char alt = snv.getNormal().getALT();
					String key = ref + "to" + alt;
					key = GenotypeKeyUtils.aggrigateKeys(key);
					SNPDepthCounter counter = null;
					if (mutationCounter.containsKey(key)) {
						counter = mutationCounter.get(key);	
						counter.reg(snv.getTumor().getTotalcnt());
						total++;
					}
					
				}
				
			}			                      
			
		}

	}

	
	public float getSNVRatio(String key){
		//6 allele
		
		SNPDepthCounter counter = null;
		key = GenotypeKeyUtils.aggrigateKeys(key);
		if (mutationCounter.containsKey(key)) {
			
			counter = mutationCounter.get(key);
			
		}
		int cnt = 1;
		if(counter!=null){
			cnt = cnt+counter.getTotal();
		}
		int total_t = total;
		if(total<cnt){
			total_t = total+6;
		}
		double r = (float)((double)cnt/(double)(total+total_t));
		double ret = r/0.1666;
		return (float)ret;
		
	}
	
	
	
}
