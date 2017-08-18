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

import java.util.Iterator;
import java.util.Map;

import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class XYSeriesExtract {

	Map<Float, DataHolderByCN> map;


	
	
	public XYSeriesExtract(Map<Float, DataHolderByCN> map) {
		this.map = map;
	}

	
	public IntervalXYDataset getXYSeriesForLogdist(float degree,
			boolean normal) {



		DataHolderByCN data = map.get(degree);
		Map<String,float[]> map = null;
		if(normal){
			map = data.getNormalLogdist();
		}else{
			map = data.getTumorLogdist();
		}
		
				
		Iterator<String> keys = map.keySet().iterator();
		XYSeriesCollection dataset = new XYSeriesCollection();
		while(keys.hasNext()){
			
			String key = keys.next();
			float[] vals = map.get(key);
			XYSeries series = new XYSeries(key);
			int n= 0;
			for(float f:vals){
				series.add(n,f);
				n++;
			}
			dataset.addSeries(series);
		}			
		
		return dataset;
	}

	public IntervalXYDataset getXYSeriesAfterFinalFilter(float degree,
			boolean normal) {

		String s1 = "_mutation";
		String s2 = normal ? "_normal" : "_tumor";
		XYSeries series1 = new XYSeries((degree*2) + s1 + s2);
		DataHolderByCN data = map.get(degree);
		if(data==null)return null;
		float[] dist = null;
		if (normal) {
			dist = data.mutationDistNFinalFilter;
		} else {
			dist = data.mutationDistTFinalFilter;

		}
		int cnt = 0;
		if (dist != null) {
			int n = 0;
			
			for (float f : dist) {
				series1.add(n, f);
				n++;
				if(f>0){
					cnt++;
				}
			}
		}
		if(series1.isEmpty()||cnt==0){
			return null;
		}
		XYSeriesCollection dataset = new XYSeriesCollection(series1);
		return dataset;
	}

	public IntervalXYDataset getXYSeries(float degree, boolean SNP,
			boolean normal, boolean filter) {

		String s1 = SNP ? "_SNP" : "_mutation";
		String s2 = normal ? "_normal" : "_tumor";

		XYSeries series1 = new XYSeries((degree*2) + s1 + s2);
		DataHolderByCN data = map.get(degree);
		if(data==null)return null;
		float[] dist = null;
		if (data != null) {

			if (SNP) {

				if (normal) {
					dist = data.snpDistN;
				} else {
					dist = data.snpDistT;
				}
			} else {

				if (filter) {
					if (normal) {
						dist = data.mutationDistNFilter;
					} else {
						dist = data.mutationDistTFilter;

					}
				} else {
					if (normal) {
						dist = data.mutationDistN;
					} else {
						dist = data.mutationDistT;

					}
				}

			}

		}
		int n = 0;
		int cnt = 0;
		if (dist != null) {
			for (float f : dist) {
				series1.add(n, f);
				n++;
				if(f>0){
					cnt++;
				}
			}
		}
		if(series1.isEmpty()||cnt==0){
			return null;
		}
		XYSeriesCollection dataset = new XYSeriesCollection(series1);
		return dataset;
	}
}
