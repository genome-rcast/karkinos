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
package jp.ac.utokyo.rcast.karkinos.graph;

import java.util.List;

import jp.ac.utokyo.karkinos.ploidy.PeakPoint;
import jp.ac.utokyo.karkinos.ploidy.TheoreticalNodes;

import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.AbstractXYZDataset;
import org.jfree.data.xy.XYZDataset;

public class MyXYZDataset extends AbstractXYZDataset {

	List<TheoreticalNodes> tnode;
	List<PeakPoint> list;
	int baseploidy;
	public MyXYZDataset(List<TheoreticalNodes> tnode, List<PeakPoint> list, int baseploidy) {
		
		this.tnode =tnode;
		this.list = list;
		this.baseploidy = baseploidy;
	}
	
	
	
	public Number getX(int i, int j) {
		
		if(i==0){
			return list.get(j).getImbalanceratio();
		}else{
			if(baseploidy==2){
				return tnode.get(i-1).getImbalance2N()[j];
			}else{
				return tnode.get(i-1).getImbalance4N()[j];
			}
		}
		
	}
	
	public Number getY(int i, int j) {

		if(i==0){
			return list.get(j).getPeakpos();
		}else{
			if(baseploidy==2){
				return tnode.get(i-1).getDistTo2N()[j];
			}else{
				return tnode.get(i-1).getDistTo4N()[j];
			}
		}
		
		
	}
	
	public Number getZ(int i, int j) {

		if(i==0){
			double mug = list.get(j).getPeakmagnitude();
			if(mug<0.01){
				mug = 0.01;
			}
			return (mug*500*0.001);
		}else{
			return (j*0.001);
		}
	}
	
	@Override
	public int getSeriesCount() {
		//
		return tnode.size()+1;
	}
	public int getItemCount(int i) {
		
		if(i==0){
			return list.size();
		}else{
			return 99;
		}		

	}
	
	@Override
	public Comparable getSeriesKey(int i) {
		
		// 
		if(i==0){
			return "observed";
			
		}else{
			return tnode.get(i-1).getID();
		}
	
	}


}
