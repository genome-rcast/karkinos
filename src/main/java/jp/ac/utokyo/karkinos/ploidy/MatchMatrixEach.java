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

import java.util.LinkedHashMap;
import java.util.Map;

public class MatchMatrixEach {


	Map<String, double[]> nodecounter;
	double sumdist;
	double sumratio;
	int numhit; 
	
	public int getNumhit() {
		return numhit;
	}
	public void setNumhit(int numhit) {
		this.numhit = numhit;
	}
	public double getSumratio() {
		return sumratio;
	}
	public void setSumratio(double sumratio) {
		this.sumratio = sumratio;
	}
	int purity;
	
	public int getPurity() {
		return purity;
	}
	public void setPurity(int purity) {
		this.purity = purity;
	}

	public Map<String, double[]> getNodecounter() {
		return nodecounter;
	}
	public void setNodecounter(Map<String, double[]> nodecounter) {
		this.nodecounter = nodecounter;
	}
	public double getSumdist() {
		return sumdist;
	}
	public void setSumdist(double sumdist) {
		this.sumdist = sumdist;
	}
	
	public String getHitnodes() {
		
		StringBuffer sb = new StringBuffer();
		if(nodecounter!=null){
			for(String s:nodecounter.keySet()){
				sb.append(",("+s+")");
			}
		}
		return sb.toString();
	}
	
	
}
