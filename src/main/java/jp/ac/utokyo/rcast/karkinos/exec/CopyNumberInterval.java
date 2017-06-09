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
package jp.ac.utokyo.rcast.karkinos.exec;

public class CopyNumberInterval implements java.io.Serializable{

	public boolean isSupportbyAllelic() {
		return supportbyAllelic;
	}

	public void setSupportbyAllelic(boolean supportbyAllelic) {
		this.supportbyAllelic = supportbyAllelic;
	}

	public void setVaridated(boolean varidated) {
		this.varidated = varidated;
	}

	public void setSnpclrrel(float snpclrrel) {
		this.snpclrrel = snpclrrel;
	}

	public void setNoSNP(int noSNP) {
		this.noSNP = noSNP;
	}
	String chr;
	int start;
	int end;
	float copynumber;
	boolean recurrent = false;
	boolean allelic = false;
	float aaf=0;
	float baf=0;
	
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

	public boolean isAllelic() {
		return allelic;
	}

	public void setAllelic(boolean allelic) {
		this.allelic = allelic;
	}

	public boolean isRecurrent() {
		return recurrent;
	}

	public void setRecurrent(boolean recurrent) {
		this.recurrent = recurrent;
	}
	boolean varidated = true;
	boolean hdelation = false;
	public boolean isHdelation() {
		return hdelation;
	}

	public void setHdelation(boolean hdelation) {
		this.hdelation = hdelation;
	}
	float snpclrrel;
	int noSNP;
	
	boolean supportbyAllelic;
	
	public String getKey(){
		return chr+":"+start+"-"+end;
	}
	
	public String getChr() {
		return chr;
	}
	public void setChr(String chr) {
		this.chr = chr;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public float getCopynumber() {
	
		return copynumber;
		
	}
	public String getRejected() {
		
		if(varidated){
			return "false";
		}else{
			if(noSNP<50){
				return "n.v";
			}
			return "true";
		}
	}
	public float getSnpclrrel() {
		return snpclrrel;
	}
	public int getNoSNP() {
		return noSNP;
	}
	public void setCopynumber(float copynumber2) {
		this.copynumber = copynumber2;
	}
	
	
	
}
