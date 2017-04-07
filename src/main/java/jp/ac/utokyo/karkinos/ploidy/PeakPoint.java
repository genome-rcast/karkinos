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

public class PeakPoint implements Comparable{
	
	float peakpos;
	float imbalanceratio;
	float peakmagnitude;
	public float getPeakpos() {
		return peakpos;
	}
	public void setPeakpos(float peakpos) {
		this.peakpos = peakpos;
	}
	public float getImbalanceratio() {
		return imbalanceratio;
	}
	public void setImbalanceratio(float imbalanceratio) {
		this.imbalanceratio = imbalanceratio;
	}
	public float getPeakmagnitude() {
		return peakmagnitude;
	}
	public void setPeakmagnitude(float peakmagnitude) {
		this.peakmagnitude = peakmagnitude;
	}
	
	public int compareTo(Object o) {
		PeakPoint p2 = (PeakPoint)o;
		return p2.peakmagnitude>peakmagnitude?1:-1;
	}

}
