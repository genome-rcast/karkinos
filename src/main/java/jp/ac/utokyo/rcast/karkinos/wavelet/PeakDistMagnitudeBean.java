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
package jp.ac.utokyo.rcast.karkinos.wavelet;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

public class PeakDistMagnitudeBean {
	

	Map<Float,Float> map = new LinkedHashMap<Float,Float>();
	
	public float getDist() {
		
		Iterator<Float> ite = map.keySet().iterator();
		double total = 0;
		double summag = 0;
		while(ite.hasNext()){
			
			float dist = ite.next();
			float mag = map.get(dist);
			total = total + (dist*mag);
			summag = summag + mag;
			
		}		
		return (float)(total/summag);
	}
	
	public Map<Float, Float> getMap() {
		return map;
	}
	public void setMap(Map<Float, Float> map) {
		this.map = map;
	}
	
	public void reg(float pd, float mg) {
		
		map.put(pd, mg);
		
	}
	
}
