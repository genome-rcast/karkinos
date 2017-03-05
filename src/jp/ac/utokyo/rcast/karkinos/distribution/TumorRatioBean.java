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


public class TumorRatioBean {

	float tumorratio;
	float observedratio;
	float sd;
	int number;
	float correl = 0;
	
	public float getCorrel() {
		return correl;
	}
	public void setCorrel(float correl) {
		this.correl = correl;
	}
	public float getTumorratio() {
		return tumorratio;
	}
	public void setTumorratio(float tumorratio) {
		this.tumorratio = tumorratio;
	}
	public float getObservedratio() {
		return observedratio;
	}
	public void setObservedratio(float observedratio) {
		this.observedratio = observedratio;
	}
	public float getSd() {
		return sd;
	}
	public void setSd(float sd) {
		this.sd = sd;
	}
	public int getNumber() {
		return number;
	}
	public void setNumber(int number) {
		this.number = number;
	}
	int mode;
	public void setMode(int i) {
		mode = i;		
	}
	public int getMode(){
		return mode;
	}
	public String getModeStr() {
		
		if(mode==0){
			return "from HetroSNP distribution";
		}else if(mode==1){
			return "from Loss/Gain base line";
		}
		return "minimal tumor ratio";
	}
	
}
