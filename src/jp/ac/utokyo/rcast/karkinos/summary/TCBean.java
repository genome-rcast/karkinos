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
package jp.ac.utokyo.rcast.karkinos.summary;

public class TCBean {

	float tumorratio = 0;
	float sd = 0;
	public float correl;
	public int nosnp;
	public int sourceflg;
	public int takefrom;
	String name;
	public float ploidy;
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public float getTumorratio() {
		return tumorratio;
	}
	public void setTumorratio(float tumorratio) {
		this.tumorratio = tumorratio;
	}
	public float getSd() {
		return sd;
	}
	public void setSd(float sd) {
		this.sd = sd;
	}
	public float getCorrel() {
		return correl;
	}
	public void setCorrel(float correl) {
		this.correl = correl;
	}
	public int getNosnp() {
		return nosnp;
	}
	public void setNosnp(int nosnp) {
		this.nosnp = nosnp;
	}
	public int getSourceflg() {
		return sourceflg;
	}
	public void setSourceflg(int sourceflg) {
		this.sourceflg = sourceflg;
	}
	public int getTakefrom() {
		return takefrom;
	}
	public void setTakefrom(int takefrom) {
		this.takefrom = takefrom;
	}
	
}
