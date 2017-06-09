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

import java.util.ArrayList;
import java.util.List;

public class SNPbean {
	
	String ref;
	String alt;
	String id;
	int cnt = 1;
	List<Float> l = new ArrayList<Float>();
	
	public SNPbean(String[] sa) {
		
		id = sa[2];
		ref = sa[3];
		alt = sa[4];
		String is = sa[7];
		float af = getAF(is);
		l.add(af);
	}
	private float getAF(String is) {
		int idxs = is.indexOf("AF=");
		String af = is.substring(idxs+3,idxs+7);
		af = af.replaceAll(",", "");
		return Float.parseFloat(af);
	}
	public String getRef() {
		return ref;
	}
	public void setRef(String ref) {
		this.ref = ref;
	}
	public String getAlt() {
		return alt;
	}
	public void setAlt(String alt) {
		this.alt = alt;
	}
	public int getCnt() {
		return cnt;
	}
	public void setCnt(int cnt) {
		this.cnt = cnt;
	}
	public float getMean(){
		
		float sum = 0;
		for(Float f:l){
			sum = sum + f;
		}
		float mean = sum/l.size();
		return mean;
	}
	
	public void inc(String[] sa){
		cnt++;
		String is = sa[7];
		float af = getAF(is);
		l.add(af);
	}

}
