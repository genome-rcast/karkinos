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

import java.util.Map.Entry;
import java.util.Set;

import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class XYSeriesForReadsDepth {
	
	ReadsSummary readSummary;
	public XYSeriesForReadsDepth(ReadsSummary readSummary){
		
		//
		this.readSummary = readSummary;
		
		
	}

	public DefaultCategoryDataset getXYSeries(){

		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		add(dataset,readSummary.getNormalDepth(),"normal");
		add(dataset,readSummary.getTumorDepth(),"tumor");
		return dataset;
	}
	
	private void add(DefaultCategoryDataset dataset, DepthCounter depth,String label) {
		
		Set<Entry<Integer,CounterA>> set=depth.getMap().entrySet();
		for(Entry<Integer,CounterA> e:set){
			//
			String s = "";
			int n = e.getKey();
			if(n>=1000){
				s= "1000 -";
//			}else if(n>=300){
//				s=""+n+"-"+(n+99);
			}else if(n>=10){
				s=""+n+"-"+(n+9);
			}else if(n>=1){
				s=""+n+"-"+(n+8);
			}else{
				s=""+n;
			}
			dataset.addValue(e.getValue().cnt,label,s);
		}		
		
	}

	public DefaultCategoryDataset getTumorXYSeries(){
		
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		add(dataset,readSummary.getTumorDepth(),"tumor");
		return dataset;
	}


	
}
