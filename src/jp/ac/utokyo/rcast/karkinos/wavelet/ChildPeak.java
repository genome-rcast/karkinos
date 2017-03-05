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

import java.util.HashSet;
import java.util.Set;

public class ChildPeak {
	
	public Set<String> getChrom() {
		return chrom;
	}
	double u;
	double v;
	double r;
	double peakdist;
	Set<String> chrom = new HashSet<String>();
	
	float aaf;
	float baf;
	float taf;
	
	public float getAaf() {
		return aaf;
	}
	public void setAaf(float aaf) {
		this.aaf = aaf;
	}
	public float getBaf() {
		return baf;
	}
	public void setBaf(float baf) {
		this.baf = baf;
	}
	public float getTaf() {
		return taf;
	}
	public void setTaf(float taf) {
		this.taf = taf;
	}
	public double getPeakdist() {
		return peakdist;
	}
	public void setPeakdist(double peakdist) {
		this.peakdist = peakdist;
	}
	public double getU() {
		return u;
	}
	public void setU(double u) {
		this.u = u;
	}
	public double getV() {
		return v;
	}
	public void setV(double v) {
		this.v = v;
	}
	public double getR() {
		return r;
	}
	public void setR(double r) {
		this.r = r;
	}
	public boolean isNotSexChrom() {
		
		if(chrom==null){
			return true;
		}else{
			for(String key:chrom){
				key = key.toUpperCase();
				if(key.contains("X")){
					return  false;
				}
				if(key.contains("Y")){
					return  false;
				}
			}
		}
		return true;
	}

}
