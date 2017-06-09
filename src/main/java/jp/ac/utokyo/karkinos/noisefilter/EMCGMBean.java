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
package jp.ac.utokyo.karkinos.noisefilter;

public class EMCGMBean {
	
	double like0;
	double like1;
	double p0;
	double p1;
	double value;

	
	public double getValue() {
		return value;
	}
	public void setValue(double value) {
		this.value = value;
	}
	public double getLike0() {
		return like0;
	}
	public void setLike0(double like0) {
		this.like0 = like0;
	}
	public double getLike1() {
		return like1;
	}
	public void setLike1(double like1) {
		this.like1 = like1;
	}
	public double getP0() {
		return p0;
	}
	public void setP0(double p0) {
		this.p0 = p0;
	}
	public double getP1() {
		return p1;
	}
	public void setP1(double p1) {
		this.p1 = p1;
	}


}
