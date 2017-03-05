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
package jp.ac.utokyo.rcast.karkinos.summary.cnvbygene;

public class CNVEachBean {

	
	String chr;
	int pos;
	float highAlleleAdj;
	float lowAlleleAdj;
	float highSmooth;
	float lowSmooth;
	float highHMM;
	float lowHMM;
	
	int start;
	int end;
	
	float ratioAdj;
	float smoothedRatio;
	float AalleleCopyNumber;
	float BalleleCopyNumber;
	float copyNumber;
	
	public void setCopyNumber(float copyNumber) {
		this.copyNumber = copyNumber;
	}

	public CNVEachBean(String[] sa, boolean total) {
			
		
		
		String chr0 = sa[0];
		if (!chr0.contains("chr")) {
			chr0 = "chr" + chr0;
		}
		chr = chr0;
		pos = Integer.parseInt(sa[1]);
		
		if(total){
			
			start = pos;
			end = Integer.parseInt(sa[2]);
			
			ratioAdj= Float.parseFloat(sa[6]);
			smoothedRatio= Float.parseFloat(sa[7]);
			AalleleCopyNumber= Float.parseFloat(sa[8]);
			BalleleCopyNumber= Float.parseFloat(sa[9]);
			copyNumber= Float.parseFloat(sa[10]);
			
		}else{
			

			highAlleleAdj= Float.parseFloat(sa[5]);
			lowAlleleAdj= Float.parseFloat(sa[6]);
			highSmooth = Float.parseFloat(sa[7]);
			lowSmooth = Float.parseFloat(sa[8]);
			highHMM = Float.parseFloat(sa[9]);
			lowHMM = Float.parseFloat(sa[10]);
			
		}
		
	}

	public float getCopyNumber() {
		return copyNumber;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public float getHighAlleleAdj() {
		return highAlleleAdj;
	}

	public float getLowAlleleAdj() {
		return lowAlleleAdj;
	}

	public float getHighSmooth() {
		return highSmooth;
	}

	public float getLowSmooth() {
		return lowSmooth;
	}

	public float getHighHMM() {
		return highHMM;
	}

	public float getLowHMM() {
		return lowHMM;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public float getRatioAdj() {
		return ratioAdj;
	}

	public float getSmoothedRatio() {
		return smoothedRatio;
	}

	public float getAalleleCopyNumber() {
		return AalleleCopyNumber;
	}

	public float getBalleleCopyNumber() {
		return BalleleCopyNumber;
	}

}
