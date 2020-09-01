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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;

public class DepthCounter implements java.io.Serializable {
	Map<Integer, CounterA> map = new TreeMap<Integer, CounterA>();
	List<Interval> lowcoverageList = new ArrayList<Interval>();

	double sumdepth = 0;
	long total = 0;

	double sumcdsdepth = 0;
	long totalcds;

	double sumontagdepth = 0;
	long totalontag=0;

	long over20x;
	long over10x;

	public float getMeanOntagDepth() {
		double mean = sumontagdepth / (double) totalontag;
		return (float) mean;
	}

	public float getMeanDepth() {
		double mean = sumdepth / (double) total;
		return (float) mean;
	}

	public float getMeanCDSDepth() {
		double mean =  sumcdsdepth  / (double) totalcds;
		return (float) mean;
	}

	public float getOver10X() {
		double mean =  over10x  / (double) totalcds;
		return (float) (mean*100);
	}

	public float getOver20X() {
		double mean =  over20x  / (double) totalcds;
		return (float) (mean*100);
	}

	public void add(String chr, int pos, int depth, int ontarget) {
//		if(ontarget==1){
//		 System.out.println("ontag=" + ontarget +" " + chr +" " + pos  +" " + depth);
//		}
		if (depth < KarkinosProp.mindepth) {
			boolean extendInterval = false;
			if (lowcoverageList.size() > 0) {
				int lastidx = lowcoverageList.size() - 1;
				Interval lastiv = lowcoverageList.get(lastidx);
				extendInterval = lastiv.extendInterval(chr, pos, depth);
				if (!extendInterval) {
					Interval iv = new Interval(chr, pos, depth);
					lowcoverageList.add(iv);
				}
			} else {
				Interval iv = new Interval(chr, pos, depth);
				lowcoverageList.add(iv);
			}
		}
		sumdepth = sumdepth + depth;
		total++;
		int depthkey =0;
		if (depth >= 1000) {
			depthkey = 1000;
			// }else if(depth>=300){
			// depth = (depth/100)*100;
		} else if (depth >= 10) {
			depthkey = (depth / 10) * 10;
		} else if (depth > 0) {
			depthkey = 1;
		}

		CounterA counter = null;
		if (!map.containsKey(depthkey)) {
			counter = new CounterA();
			map.put(depthkey, counter);
		} else {
			counter = map.get(depthkey);
		}
		counter.inc();
		
		if(ontarget==1){
			sumontagdepth = sumontagdepth +depth;
			totalontag++;
		}

		if(onCDS(chr,pos)&&ontarget==1){
			sumcdsdepth = sumcdsdepth + depth;
			 totalcds++;
			 if(depth>=10){
				 over10x++;
			 }
			 if(depth>=20){
				 over20x++;
			 }
		}
	}

	private boolean onCDS(String chr, int pos) {
		if(ge!=null){
			return ge.onCDS(chr,pos);
		}
		return false;
	}

	public List<Interval> getLowcoverageList() {
		return lowcoverageList;
	}

	public long getTotal() {
		return total;
	}

	public Map<Integer, CounterA> getMap() {
		return map;
	}

	public void merge(DepthCounter dc) {
		sumdepth = sumdepth + dc.sumdepth;
		total = total + dc.total;

		sumcdsdepth = sumcdsdepth +dc.sumcdsdepth;
		totalcds = totalcds + dc.totalcds;

		sumontagdepth = sumontagdepth +dc.sumontagdepth;
		totalontag = totalontag + dc.totalontag;

		over20x  = over20x +dc.over20x;
		over10x  = over10x + dc.over10x;

		try {
			Set<Entry<Integer, CounterA>> set = dc.getMap().entrySet();
			for (Entry<Integer, CounterA> et : set) {
				int key = et.getKey();
				CounterA value = et.getValue();
				if (map.containsKey(key)) {
					CounterA thisval = map.get(key);
					thisval.cnt = thisval.cnt + value.cnt;
				} else {
					map.put(key, value);
				}
			}
		} catch (Exception ex) {
		}
		try{
			lowcoverageList.addAll(dc.getLowcoverageList());
		}catch(Exception ex){
		}
	}

	GeneExons ge;

	public void setGeneExons(GeneExons ge) {
		this.ge = ge;	
	}
}
