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
package jp.ac.utokyo.rcast.karkinos.annotation;

public class DbSNPBean implements java.io.Serializable{

	String[] sa;
	int mode;
	float freq =0;
	boolean valid = true;
	boolean cosmic;
	boolean cosmicvalid;
	int cosmiccount;
	
	public boolean isCosmicvalid() {
		return cosmicvalid || cosmiccount>=8;
	}
	public void setCosmicvalid(boolean cosmicvalid) {
		this.cosmicvalid = cosmicvalid;
	}
	public int getCosmiccount() {
		return cosmiccount;
	}
	public void setCosmiccount(int cosmiccount) {
		this.cosmiccount = cosmiccount;
	}

	int cnt = 1;
	String varidationStr;
	
	public String getVaridationStr() {
		return varidationStr;
	}
	public void setVaridationStr(String varidationStr) {
		this.varidationStr = varidationStr;
	}
	public boolean isCosmic() {
		return cosmic;
	}
	public boolean isCosmicHigh() {
		return isCosmic()&& isValid() && cnt > 1;
	}
	

	public void setCosmic(boolean cosmic) {
		this.cosmic = cosmic;
	}

	public boolean isValid() {
		return valid;
	}

	public void setValid(boolean valid) {
		this.valid = valid;
	}

	public float getFreq() {
		return freq;
	}

	public void setFreq(float freq) {
		this.freq = freq;
	}

	public int getMode() {
		return mode;
	}

	public void setMode(int mode) {
		this.mode = mode;
	}

	public void setData(String[] _sa) {
		sa = _sa;
	}

	public String toStr() {

		StringBuffer sb = new StringBuffer();
		for (String s : sa) {
			sb.append(s+"\t");
		}
		return sb.toString();
	}

	public String getInfo() {
		if(mode==0){
			String id = sa[4];
			if(id.length()==1){
				id =".";
			}
			return id;			
		}else if(mode==DbSNPAnnotation.MODEcosmic){	
			return "cosmic_"+sa[3];
		}else if(mode==DbSNPAnnotation.MODE1000g){	
			return "1000g_"+sa[5];
		}else{
			return "exomeSNP "+ freq;
		}
	}

	public void inc() {
		cnt++;		
	}
	
}
