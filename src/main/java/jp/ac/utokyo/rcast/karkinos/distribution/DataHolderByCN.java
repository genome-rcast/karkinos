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
package jp.ac.utokyo.rcast.karkinos.distribution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.filter.FilterResult;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;
import jp.ac.utokyo.rcast.karkinos.utils.GenotypeKeyUtils;

public class DataHolderByCN implements java.io.Serializable {

	List<SNVHolder> hetroSNPlist = new ArrayList<SNVHolder>();
	float[] snpDistN = new float[100];
	float[] snpDistT = new float[100];

	List<SNVHolder> mutationSNPlist = new ArrayList<SNVHolder>();
	float[] mutationDistN = new float[100];
	float[] mutationDistT = new float[100];
	public float[] mutationDistNFilter = new float[100];
	public float[] mutationDistTFilter = new float[100];
	public float[] mutationDistNFinalFilter = new float[100];
	public float[] mutationDistTFinalFilter = new float[100];
	public float[] mutationDistTFinalFilter100xdepth = new float[100];

	
	SummaryStatistics tumorhetrosnpMeanSd = new SummaryStatistics();
			
	public void add(SNVHolder holder) {

		//
		float normalr = holder.getNormal().getRatio();
		float tumorr = holder.getTumor().getRatio();
		
		
		if (holder.isHetroSNP()) {

			hetroSNPlist.add(holder);
			addup(snpDistN, normalr);
			if(tumorr > 0.02 && tumorr < 0.98){
				addup(snpDistT, tumorr);
			}
			float r = tumorr * 100;
			if (r > 25 && r < 75) {
				tumorhetrosnpMeanSd.addValue(r);				
			}
			if(!holder.getTumor().isIndel()){
				addLogdist(holder,false);
			}
			if(!holder.getNormal().isIndel()){
				addLogdist(holder,true);
			}

		} else if ((holder.getFlg() == PileUP.SomaticMutation)
				|| (holder.getFlg() == PileUP.TumorINDEL)) {

			if (holder.getDbSNPbean() != null) {
				return;
			}

			if (tumorr > 0) {
				mutationSNPlist.add(holder);
				addup(mutationDistN, normalr);
				addup(mutationDistT, tumorr);

				if (holder.getFilterResult() != null
						&& holder.getFilterResult().isPassFilter()) {
					addup(mutationDistNFilter, normalr);
					addup(mutationDistTFilter, tumorr);
					if (CalcUtils.pass2(holder.getFilterResult())) {
						addup(mutationDistNFinalFilter, normalr);
						addup(mutationDistTFinalFilter, tumorr);
						if(holder.getTumor().getTotalcnt()>100){
							addup(mutationDistTFinalFilter100xdepth, tumorr);
						}
					}

				}

			}
		}		
		
		
	}

	
	public Map<String, float[]> getNormalLogdist() {
		return normallogdist;
	}
	public Map<String, float[]> getTumorLogdist() {
		return tumorlogdist;
	}
	Map<String,float[]> normallogdist = new TreeMap<String,float[]>();
	Map<String,float[]> tumorlogdist = new TreeMap<String,float[]>();
	
	private void addLogdist(SNVHolder holder,boolean normal) {
		////
		char ref = holder.getNormal().getGenomeR();
		char alt = holder.getNormal().getALT();
		String key = ref + "to" + alt;
		key = GenotypeKeyUtils.aggrigateKeys(key);
		if(key.equals(""))return;
		////
		double mutationlog = 0;
		if(normal){
			mutationlog = holder.getNormal().getMutationLogLikeHood()+50;
		}else{
			mutationlog = holder.getTumor().getMutationLogLikeHood()+50;
		}
		Map<String,float[]> map = normal?normallogdist:tumorlogdist;
		////
		float[] ary = null;
		if(map.containsKey(key)){
			
			//
			ary = map.get(key);
			
		}else{
			
			ary = new float[100];
			map.put(key, ary);
			
		}
		//
		int idx = getIndex(mutationlog);
		ary[idx] = ary[idx]+1;
				
	}
	
	private int getIndex(double d) {
		//d = d*100;
		int n = (int)d;
		if(n<0){
			n=0;
		}
		if(n>99){
			n=99;
		}
		return n;
	}

	public SummaryStatistics getTumorhetrosnpMeanSd() {
		return tumorhetrosnpMeanSd;
	}

	private void addup(float[] dist, float ratio) {

		//
		int idx = Math.round(ratio * 100);
		if (idx > 99)
			return;
		dist[idx] = (dist[idx] + 1);

	}

}
