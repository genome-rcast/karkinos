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

import java.util.HashMap;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class MutationCounter implements java.io.Serializable{
	
	int total=0;
	//
	int lowdepthcnt=0;
	int x2count=0;
	int x3count=0;
	int x4count=0;
	int highdepthcnt=0;
	
	int lowdepthcntH=0;
	int x2countH=0;
	int x3countH=0;
	int x4countH=0;
	int highdepthcntH=0;
	
	public void regStat(int cn,float ratio,float oddsratio,int depth){
		
		//
		total++;
		SummaryStatsHolder ss = null;
		if(!data.containsKey(cn)){
			ss = new SummaryStatsHolder();
			data.put(cn,ss);			
		}else{
			ss=data.get(cn);
		}
		//
		ss.addValue(ratio,oddsratio);
		//
		if(depth > (KarkinosProp.mindepth*5)){
			highdepthcnt++;	
		}else if(depth > (KarkinosProp.mindepth*4)){
			x4count++;	
		}else if(depth > (KarkinosProp.mindepth*3)){
			x3count++;	
		}else if(depth > (KarkinosProp.mindepth*2)){
			x2count++;		
		}else{
			lowdepthcnt++;
		}
		
	}
	
	public void regStat2( int depth) {
		if(depth > (KarkinosProp.mindepth*5)){
			highdepthcntH++;	
		}else if(depth > (KarkinosProp.mindepth*4)){
			x4countH++;	
		}else if(depth > (KarkinosProp.mindepth*3)){
			x3countH++;	
		}else if(depth > (KarkinosProp.mindepth*2)){
			x2countH++;		
		}else{
			lowdepthcntH++;
		}
		
	}
	
	public int getLowdepthcntH() {
		return lowdepthcntH;
	}

	public int getX2countH() {
		return x2countH;
	}

	public int getX3countH() {
		return x3countH;
	}

	public int getX4countH() {
		return x4countH;
	}

	public int getHighdepthcntH() {
		return highdepthcntH;
	}

	public int getX2count() {
		return x2count;
	}

	public int getX3count() {
		return x3count;
	}

	public int getX4count() {
		return x4count;
	}

	public int getLowdepthcnt() {
		return lowdepthcnt;
	}

	public int getHighdepthcnt() {
		return highdepthcnt;
	}

	public SummaryStatsHolder getSummaryStatsHolder(int key){
		if(data.containsKey(key)){
			
			return data.get(key);
		}
		return null;
	}
	
	public int getTotal() {
		return total;
	}

	Map<Integer,SummaryStatsHolder> data 
		= new HashMap<Integer,SummaryStatsHolder>();


	
	

}
