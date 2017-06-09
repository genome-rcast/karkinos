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
package jp.ac.utokyo.karkinos.ploidy;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.alleliccnv.SNVHolderPlusACnv;

public class ChrAllelicPeak {

	List<SNVHolderPlusACnv> list = new ArrayList<SNVHolderPlusACnv>();
	SummaryStatistics higha = new SummaryStatistics();
	SummaryStatistics lowa = new SummaryStatistics();
	AFCounter tumorSNP = new AFCounter();
	boolean even = false;
	
	public boolean isEven() {
		return even;
	}

	public void setEven(boolean even) {
		this.even = even;
	}

	public void add(SNVHolderPlusACnv snva){
		
		//
		list.add(snva);
		higha.addValue(snva.getHighera().getWtval());
		lowa.addValue(snva.getLowera().getWtval());
		//
		tumorSNP.setValue(snva.getSnv().getTumor().getRatio());
		
	}

	public AFCounter getTumorSNP() {
		return tumorSNP;
	}
	
	
	
}
