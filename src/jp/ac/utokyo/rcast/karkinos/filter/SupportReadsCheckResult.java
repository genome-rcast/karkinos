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
package jp.ac.utokyo.rcast.karkinos.filter;

import java.util.Set;

public class SupportReadsCheckResult {

	Set<Integer> filter;
	float pval4directionCheck = 0;
	float logtAdjusted=0;
	
	float ntQualityDiff=0;
	
	public float getNtQualityDiff() {
		return ntQualityDiff;
	}
	public void setNtQualityDiff(float ntQualityDiff) {
		this.ntQualityDiff = ntQualityDiff;
	}
	public float getSupportreadsBAlleleFeqquency() {
		return supportreadsBAlleleFeqquency;
	}
	public float getRefreadsBAlleleFeqquency() {
		return refreadsBAlleleFeqquency;
	}
	float supportreadsBAlleleFeqquency = 0;
	float refreadsBAlleleFeqquency = 0;
	//			
		
	public float getLogtAdjusted() {
		return logtAdjusted;
	}
	public void setSupportreadsBAlleleFeqquency(float supportreadsBAlleleFeqquency) {
		this.supportreadsBAlleleFeqquency = supportreadsBAlleleFeqquency;
	}
	public void setRefreadsBAlleleFeqquency(float refreadsBAlleleFeqquency) {
		this.refreadsBAlleleFeqquency = refreadsBAlleleFeqquency;
	}
	public void setLogtAdjusted(float logt) {
		this.logtAdjusted = logt;
	}
	public Set<Integer> getFilter() {
		return filter;
	}
	public void setFilter(Set<Integer> filter) {
		this.filter = filter;
	}
	public float getPval4directionCheck() {
		return pval4directionCheck;
	}
	public void setPval4directionCheck(float pval4directionCheck) {
		this.pval4directionCheck = pval4directionCheck;
	}
	
}